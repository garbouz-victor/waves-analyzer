"""STEP3A4 role-separated, append-only scientific stages. Historical inputs read-only."""
import hashlib
import csv
import json
import os
from pathlib import Path
import tarfile
import time
import numpy as np
from ..ch_validation_policy import (ValidationPolicy, AnalysisBudget, LINEAR_LIMITS,
    authorizations, linear_gates, run_sparse_analysis)
from ..linearized_ch import json_hash, array_hash
from ..performance import Performance, petsc_events
from . import spectral_ch as historical

RESULTS = Path("validation_results/step3a4")
POLICY = ValidationPolicy()
ACTIVE_OUTPUT = None


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def freeze_policy():
    folder = RESULTS/"policy"
    folder.mkdir(parents=True, exist_ok=True)
    record = {"policy": POLICY.record(), "policy_sha256": POLICY.fingerprint,
        "design_sha256": hashlib.sha256(Path("STEP3A4_DESIGN.md").read_bytes()).hexdigest()}
    path = folder/"policy.json"
    if path.exists():
        if read_json(path) != record:
            raise ValueError("STOP: frozen STEP3A4 policy or design changed")
    else:
        write_json(path, record)
    return record


def reserve(name):
    global ACTIVE_OUTPUT
    output = RESULTS/name
    if RESULTS.resolve() not in output.resolve().parents:
        raise ValueError("Historical/external output paths forbidden; require STEP3A4 subdirectory")
    policy = freeze_policy()
    output.mkdir(parents=True, exist_ok=False)
    ACTIVE_OUTPUT = output
    write_json(output/"status.json", {"status": "running", "stage": name})
    source = Path(__file__).resolve().parents[1]
    with tarfile.open(output/"multiphase_source.tar.gz", "x:gz") as archive:
        for path in sorted(source.rglob("*.py")):
            archive.add(path, arcname=str(path.relative_to(source)))
    from ..provenance import source_hash
    old = read_json(historical.RESULTS/"linear_operator/status.json")["provenance"]
    provenance = {"git_commit": os.environ.get("STEP3_GIT_COMMIT", "unavailable"),
        "docker_image_digest": os.environ.get("STEP3_DOCKER_IMAGE_DIGEST", "unavailable"),
        "multiphase_source_sha256": source_hash(),
        "source_archive_sha256": hashlib.sha256((output/"multiphase_source.tar.gz").read_bytes()).hexdigest(),
        "historical_input_provenance": old, **policy}
    for key in ("equilibrium_fingerprint", "perturbation_coefficient_sha256", "initial_phi_sha256",
                "mesh_fingerprint", "quadrature_degree", "mobility", "target_mass"):
        provenance[key] = old[key]
    write_json(output/"provenance.json", provenance)
    return output


def analysis_used():
    total = 0.
    for path in RESULTS.rglob("timing.json"):
        record = read_json(path)
        if record.get("analysis_budget_charged", False):
            total += record["wall_s"]
    return total


def load_reference():
    from ..ch_timestep_planner import ReducedCHReference
    op, _, operator_record, reference_record, basis, _ = historical.qualified_reference()
    V = basis["V"]
    chemical = op.mass_solve(op.H @ V)
    ref = ReducedCHReference(basis["projected_L"], V.T @ (op.H @ V),
        op.mobility*chemical.T @ (op.K @ chemical),
        float(basis["beta"])*np.eye(V.shape[1])[:, 0], V.T @ (op.M0 @ V))
    fps = {"operator": op.fingerprint,
        "equilibrium": operator_record["provenance"]["equilibrium_fingerprint"],
        "perturbation": operator_record["provenance"]["perturbation_coefficient_sha256"],
        "Krylov": reference_record["Krylov_sha256"]}
    return op, ref, V, fps, reference_record


def profile_initial():
    """No time advance: never scientific temporal evidence."""
    from dolfinx import fem
    from petsc4py import PETSc
    from ..diagnostics import Diagnostics, interface_resolution
    from ..be_work_diagnostics import BEWorkDiagnostics
    from .contact_angle import observe
    output = reserve("performance/initial_profile")
    perf, budget = Performance(), AnalysisBudget(analysis_used())
    PETSc.Log.begin()
    with perf.measure("initialization"):
        s, phi_eq, psi, initial, provenance = historical.scientific_input()
        field = fem.Function(phi_eq.function_space)
        field.x.array[:] = initial; field.x.scatter_forward()
        s.initialize(field)
        if array_hash(s.state.sub(2).collapse().x.array) != provenance["initial_phi_sha256"]:
            raise ValueError("Initial phi hash changed")
        diagnostic, work = Diagnostics(s), BEWorkDiagnostics(s)
    budget.check()
    samples = []
    for _ in range(3):
        with perf.measure("diagnostics"):
            row = diagnostic.measure()
        with perf.measure("interface_resolution"):
            resolution = interface_resolution(s)
        with perf.measure("contact_geometry"):
            geometry = observe(s)  # includes its own resolution; disclosed in result
        with perf.measure("BE_work"):
            audit = work.measure(0., row, row)
        sample = {"row": row, "resolution": resolution, "geometry": geometry, "work": audit}
        with perf.measure("filesystem_write"):
            with (output/"samples.jsonl").open("a") as stream:
                stream.write(json.dumps(sample)+"\n")
        samples.append(sample)
        budget.check()
    timing = {**perf.snapshot(), "analysis_budget_charged": True,
        "PETSc": petsc_events(), "time_advance": False,
        "geometry_timer_includes_duplicate_resolution": True}
    write_json(output/"timing.json", timing)
    result = {"status": "complete", "qualification_status": "profile_only",
        "provenance": provenance, "initial": samples[0], "timing": timing}
    write_json(output/"status.json", result)
    print(json.dumps({"stage": "initial_profile", "timing": timing["categories"]}), flush=True)
    return result


def sparse_stage(optimized=False):
    from petsc4py import PETSc
    from ..ch_timestep_planner import validate_plan, verify_full_sparse_schedule
    folder = "optimized_plan" if optimized else "baseline_plan"
    output = reserve("linear_sparse/"+folder)
    perf, budget = Performance(), AnalysisBudget(analysis_used())
    PETSc.Log.begin()
    try:
        with perf.measure("initialization"):
            op, ref, V, fps, ref_record = load_reference()
            path = (RESULTS/"schedule_optimization/optimized_plan.json" if optimized else
                    historical.RESULTS/"stiff_bump/planner/timestep_plan.json")
            plan = read_json(path)
            validate_plan(plan, fps)
            reduced_ok = all(linear_gates(plan["predicted_errors"]).values())
            op.resolvent_backend, op.performance = "petsc_mumps", perf
        write_json(output/"input_plan.json", plan)
        steps = sum(b["steps"] for b in plan["blocks"])
        def monitor(row):
            with perf.measure("filesystem_write"):
                with (output/"progress.jsonl").open("a") as stream:
                    stream.write(json.dumps(row)+"\n")
                write_json(output/"timing_progress.json", perf.snapshot())
            print("Full sparse "+folder+" "+json.dumps(row), flush=True)
            # A bounded measured forecast, not the historical nonlinear timing.
            records = perf.records
            done = len(records["linear_solve"])
            per_step = sum(r["wall_s"] for key in ("linear_solve", "diagnostics")
                           for r in records[key])/max(done, 1)
            factors = records["factorization_setup"]
            per_factor = sum(r["wall_s"] for r in factors)/max(len(factors), 1)
            budget.check(per_step*(steps-done)+per_factor*(len(plan["blocks"])-len(factors)))
        verification = run_sparse_analysis(lambda: verify_full_sparse_schedule(
            op, ref, V, plan["blocks"], monitor, performance=perf, check_budget=budget.check),
            reduced_pass=reduced_ok, analysis_budget=budget)
        gates = linear_gates(verification, full_sparse=True)
        with perf.measure("filesystem_write"):
            write_json(output/"verification.json", verification)
        status = {"status": "complete", "qualification_status": "passed" if all(gates.values()) else "failed",
            "benchmark": "full sparse LINEAR BE, not nonlinear isolated CH or CHNS",
            "steps": steps, "gates": gates, "plan_sha256": plan["plan_sha256"],
            "input_fingerprints": fps, "rational_reference_record": ref_record,
            "initial_D2": op.power(V @ ref.seed), "initial_E2": op.energy(V @ ref.seed),
            "errors": {k: v for k, v in verification.items() if k not in ("history", "blocks")},
            "reduced_errors": plan["predicted_errors"],
            "authorizations": authorizations(reduced=reduced_ok, sparse=all(gates.values()))}
        write_json(output/"status.json", status)
        print(json.dumps({"stage": folder, "gates": gates, "errors": status["errors"]}), flush=True)
        return status
    except Exception as exc:
        write_json(output/"status.json", {"status": "stopped", "qualification_status": "failed",
            "reason": str(exc), "authorizations": authorizations()})
        raise
    finally:
        timing = {**perf.snapshot(), "analysis_budget_charged": True, "PETSc": petsc_events(),
            "distinct_factors": len({r["dt"] for r in perf.records["factorization_setup"]})}
        write_json(output/"timing.json", timing)


def optimizer_stage():
    from ..ch_schedule_optimizer import optimize, original_from_uniform
    from ..ch_timestep_planner import validate_plan
    baseline_result = read_json(RESULTS/"linear_sparse/baseline_plan/status.json")
    if baseline_result["qualification_status"] != "passed":
        raise ValueError("STOP: baseline full-sparse prerequisite failed")
    output = reserve("schedule_optimization")
    perf, budget = Performance(), AnalysisBudget(analysis_used())
    optimizer_cpu_start = time.process_time()
    candidates = output/"candidates"; candidates.mkdir()
    try:
        with perf.measure("initialization"):
            op, ref, V, fps, ref_record = load_reference()
            baseline = read_json(historical.RESULTS/"stiff_bump/planner/timestep_plan.json")
            validate_plan(baseline, fps)
            original = original_from_uniform(baseline)
        def candidate(row):
            budget.check()
            if time.process_time()-optimizer_cpu_start > POLICY.optimizer_cpu_s:
                raise RuntimeError("STOP: total optimizer CPU guard including reproducibility repeat")
            with perf.measure("filesystem_write"):
                write_json(candidates/f"candidate_{row['candidate']:04d}.json", row)
            if row["action"].startswith("accept") or row["action"] == "original":
                print("Optimizer "+json.dumps(row), flush=True)
        with perf.measure("reduced_optimization"):
            found = optimize(ref, original, fps, reference_qualified=True, candidate_callback=candidate)
        with perf.measure("reproducibility_repeat"):
            def repeat_guard(row):
                budget.check()
                if time.process_time()-optimizer_cpu_start > POLICY.optimizer_cpu_s:
                    raise RuntimeError("STOP: total optimizer CPU guard including reproducibility repeat")
            repeated = optimize(ref, original, fps, reference_qualified=True, candidate_callback=repeat_guard)
        reproducible = found["plan_sha256"] == repeated["plan_sha256"]
        if not reproducible:
            raise ValueError("STOP: optimizer semantic SHA not reproducible")
        write_json(output/"optimization.json", found)
        use_found = found["qualification_status"] == "linear_reduced_pass" and found["total_steps"] < baseline["total_steps"]
        # Provenance/timing of optimization is separate from its deterministic history.
        plan = {k: v for k, v in (found if use_found else baseline).items() if k not in (
            "plan_sha256", "nonlinear_execution_authorized", "cost_gate", "cost_estimate",
            "full_sparse_verification", "full_sparse_verification_required")}
        plan.update(authorizations=authorizations(reduced=True),
            plan_selection="optimized" if use_found else "full-sparse-verified baseline fallback",
            baseline_plan_sha256=baseline["plan_sha256"],
            optimizer_reproducible=reproducible, optimizer_result_sha256=found["plan_sha256"],
            provenance=read_json(output/"provenance.json"))
        plan["plan_sha256"] = json_hash(plan)
        write_json(output/"optimized_plan.json", plan)
        budget.check()
        result = {"status": "complete", "qualification_status": "passed",
            "total_steps": plan["total_steps"], "saved_intervals": baseline["total_steps"]-plan["total_steps"],
            "reduction_percent": 100*(1-plan["total_steps"]/baseline["total_steps"]),
            "reproducible": reproducible, "plan_sha256": plan["plan_sha256"],
            "predicted_errors": plan["predicted_errors"], "requires_full_sparse": True}
        write_json(output/"status.json", result)
        print(json.dumps(result), flush=True)
        return result
    finally:
        write_json(output/"timing.json", {**perf.snapshot(), "analysis_budget_charged": True})


def isolated_stage(level=0, pilot=False):
    from dolfinx import fem
    from petsc4py import PETSc
    from ..isolated_ch import IsolatedCH
    from ..diagnostics import Diagnostics, interface_resolution
    from ..be_work_diagnostics import BEWorkDiagnostics
    from ..ch_validation_diagnostics import ValidationSchedule, WallTrace, observation_indices
    from .contact_angle import observe
    if level not in (0, 1, 2) or (pilot and level != 0):
        raise ValueError("Only level0 pilot or the three fixed refinement levels are allowed")
    plan = read_json(RESULTS/"schedule_optimization/optimized_plan.json")
    sparse = read_json(RESULTS/"linear_sparse/optimized_plan/status.json")
    if sparse.get("qualification_status") != "passed":
        raise ValueError("STOP: optimized full-sparse prerequisite failed")
    if not pilot:
        cost = read_json(RESULTS/"performance/isolated_pilot/status.json")
        if not cost.get("isolated_CH_execution_authorized", False):
            raise ValueError("STOP: isolated CH pilot has not authorized the three-level series")
        if level:
            previous_level = read_json(RESULTS/f"isolated_stiff_bump/level{level-1}/status.json")
            if (previous_level.get("qualification_status") != "level_complete" or
                    previous_level.get("optimized_plan_sha256") != plan["plan_sha256"]):
                raise ValueError("STOP: preceding matching isolated level has not completed")
    prior_isolated_wall = sum(read_json(path)["wall_s"] for path in RESULTS.rglob("timing.json")
        if "isolated_stiff_bump" in path.parts or "isolated_pilot" in path.parts)
    if prior_isolated_wall >= POLICY.isolated_series_wall_s:
        raise ValueError("STOP: total isolated series cost guard")
    name = "performance/isolated_pilot" if pilot else f"isolated_stiff_bump/level{level}"
    output = reserve(name)
    perf = Performance()
    PETSc.Log.begin()
    histories, audits, logs = [], [], []
    engine = None
    try:
        with perf.measure("initialization"):
            s, phi_eq, psi, initial, provenance = historical.scientific_input()
            eq_energy = read_json(historical.EQUILIBRIUM/"prepared.json")["free_energy"]
            field = fem.Function(phi_eq.function_space)
            field.x.array[:] = initial; field.x.scatter_forward()
            role = "isolated_CH_cost_pilot" if pilot else "isolated_CH"
            schedule = ValidationSchedule(plan, plan["input_fingerprints"], sparse, role,
                {role+"_execution_authorized": True}, factor=2**level)
            engine = IsolatedCH(s, field)
            engine.performance = perf
            if array_hash(s.state.sub(2).collapse().x.array) != provenance["initial_phi_sha256"]:
                raise ValueError("STOP: initial phi coefficient SHA mismatch")
            diagnostic, work, trace = Diagnostics(s), BEWorkDiagnostics(s), WallTrace(s)
            reference = read_json(historical.RESULTS/"krylov/rational_balanced/status.json")
            selections = observation_indices(schedule.descriptor["blocks"],
                1/reference["rho_relevant"], 4.862048987156724e-7)
            area = (s.config.x_max-s.config.x_min)*(s.config.z_max-s.config.z_min)
        nsteps = min(POLICY.isolated_pilot_steps, schedule.nsteps) if pilot else schedule.nsteps
        write_json(output/"run_plan.json", {**schedule.descriptor, "diagnostic_selection": selections,
            "pilot": pilot, "executed_intervals_target": nsteps, "role": role,
            "optimized_plan_sha256": plan["plan_sha256"]})
        write_json(output/"input_provenance.json", {**provenance,
            "optimized_plan_sha256": plan["plan_sha256"], "benchmark": "FULL NONLINEAR isolated CH; u=0 by construction"})
        def sample():
            with perf.measure("diagnostics"):
                row = diagnostic.measure()
            with perf.measure("interface_resolution"):
                resolution = interface_resolution(s)
            with perf.measure("wall_topology"):
                crossings = trace.crossings()
            row.update(dt=s.current_dt or 0., step=s.step_number,
                mass_error_to_target_domain=abs(row["phase_mass"]-provenance["target_mass"])/area,
                cells_across_transition_certified_min=resolution.get("cells_across_transition_certified_min", 0.),
                wall_crossings=len(crossings), excess_free_energy=row["E_total"]-eq_energy)
            if not np.isfinite(list(row.values())).all():
                raise RuntimeError("STOP: nonfinite isolated benchmark diagnostic")
            with perf.measure("filesystem_write"):
                with (output/"resolutions.jsonl").open("a") as stream:
                    stream.write(json.dumps(resolution)+"\n")
                with (output/"wall_topology.jsonl").open("a") as stream:
                    stream.write(json.dumps({"time": s.time, "crossings": crossings})+"\n")
            histories.append(row)
            if (row["mass_error_to_target_domain"] > POLICY.nonlinear_mass_domain or
                    not resolution["qualified_resolution"] or len(crossings) != 2):
                raise RuntimeError("STOP: isolated mass / certified resolution / contact topology gate")
            if s.step_number in selections["geometry"]:
                with perf.measure("contact_geometry"):
                    geometry = observe(s)
                if not np.allclose(crossings, sorted(c["coordinate"] for c in geometry["crossings"]), atol=1e-9, rtol=0):
                    raise RuntimeError("Cached facet topology disagrees with full contour")
                with perf.measure("filesystem_write"):
                    with (output/"geometry.jsonl").open("a") as stream:
                        stream.write(json.dumps(geometry)+"\n")
            return row
        initial_row = previous = sample()
        write_json(output/"initial.json", initial_row)
        with (output/"history.csv").open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=initial_row.keys())
            writer.writeheader(); writer.writerow(initial_row); stream.flush()
            for i in range(nsteps):
                interval = schedule.step(i)
                def snes_monitor(snes, iteration, residual):
                    with (output/"snes_iterations.jsonl").open("a") as record:
                        record.write(json.dumps({"interval": i+1, "iteration": iteration,
                            "residual": residual, "dt": interval["dt"]})+"\n")
                engine.problem.solver.cancelMonitor()
                engine.problem.solver.setMonitor(snes_monitor)
                log = engine.step(interval["dt"], interval["time"])
                logs.append(log)
                row = sample()
                if row["delta_E_interval"] > POLICY.energy_growth_J_per_m:
                    raise RuntimeError("STOP: isolated free-energy growth")
                if s.step_number in selections["BE_work"]:
                    with perf.measure("BE_work"):
                        audit = work.measure(interval["dt"], previous, row)
                    audits.append(audit)
                    with perf.measure("filesystem_write"):
                        with (output/"be_work.jsonl").open("a") as record:
                            record.write(json.dumps(audit)+"\n")
                    if (abs(audit["weak_work_defect"]) > POLICY.weak_work_J_per_m or
                            abs(audit["decomposition_roundoff"]) > POLICY.weak_work_J_per_m):
                        raise RuntimeError("STOP: isolated BE-work diagnostic gate")
                with perf.measure("filesystem_write"):
                    writer.writerow(row); stream.flush()
                    with (output/"snes.jsonl").open("a") as record:
                        record.write(json.dumps(log)+"\n")
                    if s.step_number in selections["probe_indices"] or s.step_number == nsteps:
                        np.savez_compressed(output/f"phi_n{s.step_number}.npz", phi=s.state.sub(2).collapse().x.array)
                previous = row
                if i < 3 or (i+1) % 25 == 0:
                    print(f"{role} level{level} {i+1}/{nsteps}: t={s.time:.12g} D={row['CH_dissipation']:.9g}", flush=True)
                elapsed = perf.snapshot()["wall_s"]
                if ((level == 0 and elapsed > POLICY.isolated_level0_wall_s) or
                        elapsed+prior_isolated_wall > POLICY.isolated_series_wall_s):
                    raise RuntimeError("STOP: isolated CH wall cost gate during execution")
        result = {"status": "complete", "qualification_status": "pilot_only" if pilot else "level_complete",
            "accepted_steps": nsteps, "planned_steps": schedule.nsteps, "level": level,
            "benchmark": "FULL NONLINEAR isolated CH; Ekin/Dvisc/Dslip=0 BY CONSTRUCTION",
            "optimized_plan_sha256": plan["plan_sha256"], "initial_phi_sha256": provenance["initial_phi_sha256"],
            "DeltaF": previous["E_total"]-initial_row["E_total"],
            "integrated_D_CH": previous["cumulative_CH_dissipation"],
            "final_budget_defect": previous["energy_budget_defect"],
            "max_mass_domain": max(r["mass_error_to_target_domain"] for r in histories),
            "min_certified_cells": min(r["cells_across_transition_certified_min"] for r in histories),
            "max_weak_work_defect": max((abs(a["weak_work_defect"]) for a in audits), default=0.),
            "first_BE_work": audits[0] if audits else None, "runtime_s": perf.snapshot()["wall_s"]}
        if pilot:
            timing = perf.snapshot()
            old = read_json(historical.RESULTS/"low_mode/status.json")
            old_snes = max(r["runtime_s"]/r["steps"] for r in old["temporal"])
            measured_snes = timing["categories"]["nonlinear_SNES"]["wall_s"]/nsteps
            diagnostic_mean = sum(timing["categories"].get(k, {}).get("wall_s", 0.)
                for k in ("diagnostics", "interface_resolution", "wall_topology", "filesystem_write"))/(nsteps+1)
            # Include sparse geometry/audit observations, scaled by full-series frequency.
            geometry_mean = timing["categories"]["contact_geometry"]["wall_s"]/timing["categories"]["contact_geometry"]["count"]
            audit_mean = timing["categories"]["BE_work"]["wall_s"]/timing["categories"]["BE_work"]["count"]
            initialization = timing["categories"]["initialization"]["wall_s"]
            predictions = []
            for factor in (1, 2, 4):
                from ..ch_timestep_planner import refine_blocks
                selection = observation_indices(refine_blocks(plan["blocks"], factor),
                    1/reference["rho_relevant"], 4.862048987156724e-7)
                predictions.append(initialization+(max(old_snes, measured_snes)+diagnostic_mean)*plan["total_steps"]*factor+
                    geometry_mean*len(selection["geometry"])+audit_mean*len(selection["BE_work"]))
            allowed = predictions[0] <= POLICY.isolated_level0_wall_s and sum(predictions)+result["runtime_s"] <= POLICY.isolated_series_wall_s
            result.update(measured_SNES_seconds_per_step=measured_snes,
                historical_low_mode_SNES_seconds_per_step=old_snes,
                predicted_level_wall_s=predictions, predicted_series_wall_s=sum(predictions)+result["runtime_s"],
                isolated_CH_execution_authorized=allowed)
        write_json(output/"status.json", result)
        return result
    except Exception as exc:
        record = {"status": "stopped", "qualification_status": "failed", "reason": str(exc),
            "accepted_steps": s.step_number if engine is not None else 0,
            "attempted_interval": s.step_number+1 if engine is not None else None,
            "isolated_CH_execution_authorized": False, "full_CHNS_probe_authorized": False,
            "optimized_plan_sha256": plan["plan_sha256"], "pilot": pilot}
        if engine is not None:
            snes = engine.problem.solver
            record.update(SNES_reason=snes.getConvergedReason(), SNES_residual=snes.getFunctionNorm(),
                          SNES_iterations=snes.getIterationNumber())
            np.savez_compressed(output/"failed_iterate.npz", mixed=engine.state.x.array,
                old_mixed=engine.old.x.array, accepted_phi=s.state.sub(2).collapse().x.array,
                failed_phi=engine.state.sub(0).collapse().x.array, failed_mu=engine.state.sub(1).collapse().x.array,
                initial_phi=initial, phi_eq=phi_eq.x.array)
        write_json(output/"status.json", record)
        print(json.dumps(record), flush=True)
        raise
    finally:
        write_json(output/"timing.json", {**perf.snapshot(), "analysis_budget_charged": False, "PETSc": petsc_events()})


def failed_interval_audit():
    """Read-only reassembly of a REJECTED iterate: never a retry or accepted step."""
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from ..isolated_ch import IsolatedCH
    from ..diagnostics import Diagnostics
    from ..equilibrium import _measures
    from ..linearized_ch import NonlinearCHRate
    from ..be_work_diagnostics import BEWorkDiagnostics
    pilot = RESULTS/"performance/isolated_pilot"
    failed = read_json(pilot/"status.json")
    if failed["status"] != "stopped" or failed["accepted_steps"] != 0:
        raise ValueError("This postmortem requires the stored first-interval rejection")
    output = reserve("performance/failed_interval_audit")
    perf, budget = Performance(), AnalysisBudget(analysis_used())
    with perf.measure("initialization"):
        s, phi_eq, psi, initial, provenance = historical.scientific_input()
        field = fem.Function(phi_eq.function_space)
        field.x.array[:] = initial; field.x.scatter_forward()
        iso = IsolatedCH(s, field)
        op, _, _ = historical.load_operator()
        dt = read_json(pilot/"run_plan.json")["blocks"][0]["dt"]
        iso.dt.value = dt
        diagnostic = Diagnostics(s)
        old_row = diagnostic.measure()
        with np.load(pilot/"failed_iterate.npz") as data:
            if not np.array_equal(data["initial_phi"], initial):
                raise ValueError("Failed-state input fingerprint changed")
            iso.state.x.array[:] = data["mixed"]
            iso.old.x.array[:] = data["old_mixed"]
        iso.state.x.scatter_forward(); iso.old.x.scatter_forward()
    budget.check()
    def assembled(form):
        vector = fp.assemble_vector(fem.form(form))
        result = vector.array.copy(); vector.destroy()
        return result
    with perf.measure("assembly"):
        _, pmap = iso.state.function_space.sub(0).collapse()
        _, mumap = iso.state.function_space.sub(1).collapse()
        phi, mu = ufl.split(iso.state)
        old_phi, _ = ufl.split(iso.old)
        test, _ = ufl.TestFunctions(iso.state.function_space)
        dx, _ = _measures(s)
        p = iso.state.sub(0).collapse().x.array.copy()
        pold = iso.old.sub(0).collapse().x.array.copy()
        mu_array = iso.state.sub(1).collapse().x.array.copy()
        increment = p-pold  # analysis re-expression ONLY; never assigned to any solution
        delta_field = fem.Function(phi_eq.function_space)
        delta_field.x.array[:] = increment; delta_field.x.scatter_forward()
        phase = assembled(iso.phase_form)[pmap]
        chemical = assembled(iso.chemical_form)[mumap]
        time_literal = assembled((phi-old_phi)/iso.dt*test*dx)[pmap]
        flux = assembled(s.config.mobility*ufl.inner(ufl.grad(mu), ufl.grad(test))*dx)[pmap]
        time_coefficient = op.M0 @ increment/dt
        sparse_phase = time_coefficient+s.config.mobility*(op.K @ mu_array)
        chemical_sparse = op.M0 @ mu_array-NonlinearCHRate(s, phi_eq.function_space, op).residual(p)
        stable_time_ufl = assembled(delta_field/iso.dt*ufl.TestFunction(phi_eq.function_space)*dx)
        def scalar(form):
            return float(fem.assemble_scalar(fem.form(form)))
        interpolation_difference_L2 = np.sqrt(max(0., scalar(((phi-old_phi)-delta_field)**2*dx)))
    with perf.measure("diagnostics"):
        # Publish the rejected snapshot ONLY to postmortem diagnostics, not advance().
        s.older.x.array[:] = s.state.x.array
        s.state.sub(2).interpolate(iso.state.sub(0).collapse())
        s.state.sub(3).interpolate(iso.state.sub(1).collapse())
        s.state.x.scatter_forward(); s.older.x.scatter_forward()
        s.time = dt
        rejected_row = diagnostic.measure()
        ulp = np.abs(np.spacing(pold))
        record = {"status": "complete", "purpose": "conditioning diagnosis of REJECTED stored iterate, no solve/acceptance",
            "dt": dt, "accepted_steps": 0, "source_failed_state_sha256": hashlib.sha256((pilot/"failed_iterate.npz").read_bytes()).hexdigest(),
            "stored_snapshot_note": "last evaluated Function after SNES exception; not guaranteed identical to internal last accepted Newton vector",
            "algebraic_norms_diagnostic_only": {
                "phase_literal_UFL": float(np.linalg.norm(phase)),
                "chemical_literal_UFL": float(np.linalg.norm(chemical)),
                "combined_literal_UFL": float(np.linalg.norm(np.r_[phase, chemical])),
                "phase_coefficient_difference_sparse": float(np.linalg.norm(sparse_phase)),
                "chemical_sparse": float(np.linalg.norm(chemical_sparse)),
                "literal_time_minus_coefficient_time": float(np.linalg.norm(time_literal-time_coefficient)),
                "coefficient_time_UFL_minus_sparse": float(np.linalg.norm(stable_time_ufl-time_coefficient)),
                "literal_flux_minus_sparse": float(np.linalg.norm(flux-s.config.mobility*(op.K @ mu_array))),
                "phase_form_minus_separate_time_plus_flux": float(np.linalg.norm(phase-time_literal-flux)),
                "literal_time": float(np.linalg.norm(time_literal)), "flux": float(np.linalg.norm(flux))},
            "phase_Riesz_L2_literal": op.dual_norm(phase),
            "phase_Riesz_L2_coefficient_expression": op.dual_norm(sparse_phase),
            "interpolation_subtraction_discrepancy_L2": float(interpolation_difference_L2),
            "interpolation_subtraction_discrepancy_L2_over_dt": float(interpolation_difference_L2/dt),
            "delta_phi_L2": op.norm(increment), "delta_phi_max": float(max(abs(increment))),
            "one_coefficient_ULP_over_dt_max": float(max(ulp)/dt),
            "rejected_mass_change_domain": float(abs(op.mass @ increment)/op.area),
            "no_tolerance_change": True, "no_new_time_step": True, "no_solution_correction": True,
            "SNES_atol": s.config.snes_atol, "SNES_rtol": s.config.snes_rtol,
            "initial": old_row, "rejected_state": rejected_row}
    with perf.measure("BE_work"):
        record["rejected_BE_work_diagnostic"] = BEWorkDiagnostics(s).measure(dt, old_row, rejected_row)
    record["rejected_BE_work_diagnostic"].update(accepted_interval=False,
        label="REJECTED iterate only: no nonlinear trajectory or qualification evidence")
    budget.check()
    write_json(output/"status.json", record)
    write_json(output/"timing.json", {**perf.snapshot(), "analysis_budget_charged": True})
    print(json.dumps({k: v for k, v in record.items() if k not in ("initial", "rejected_state")}), flush=True)
    return record


STAGES = {"profile": profile_initial, "baseline_sparse": sparse_stage, "optimize": optimizer_stage,
          "optimized_sparse": lambda: sparse_stage(True), "isolated_pilot": lambda: isolated_stage(pilot=True),
          "isolated_level0": lambda: isolated_stage(0), "isolated_level1": lambda: isolated_stage(1),
          "isolated_level2": lambda: isolated_stage(2), "failure_audit": failed_interval_audit}
