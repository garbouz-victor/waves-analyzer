"""Append-only STEP3A5 experiments; each command executes ONE authorized stage."""
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile
import time
import numpy as np
from ..linearized_ch import array_hash, json_hash
from ..phase_rate import VERSION, ABSOLUTE_VERSION, snes_options
from ..performance import Performance, petsc_events
from . import spectral_ch as historical

RESULTS = Path("validation_results/step3a5")
PRIOR = Path("validation_results/step3a4")
ACTIVE_OUTPUT = None
TINY_DT = 4.657657028534679e-12
MODERATE_DT = 1e-6
INITIAL_D = .11400313905309394
HASHES = {
    "initial_phi_sha256": "6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681",
    "perturbation_coefficient_sha256": "831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0",
    "equilibrium_fingerprint": "df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc",
    "optimized_plan_sha256": "5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8"}
POLICY = {"tiny_dt": TINY_DT, "moderate_dt": MODERATE_DT,
    "snes_atol": 1e-10, "snes_rtol": 1e-9, "snes_stol": 0., "snes_max_it": 30,
    "snes_type": "newtonls", "snes_linesearch_type": "bt",
    "moderate_phi_mu_D_relative": 1e-10, "moderate_energy_mass_absolute": 1e-12,
    "linear_endpoint_relative": 1e-10, "linear_increment_relative": 1e-6,
    "linear_mass_domain": 1e-11, "nonlinear_increment_relative": 1e-3,
    "mass_domain": 1e-10, "certified_cells": 8., "energy_growth_J_per_m": 1e-9,
    "weak_work_J_per_m": 1e-12, "root_cause_residual_ratio": 10.,
    "initial_guess": "unprojected semidiscrete rate; unchanged consistent mu",
    "pilot_guess": "previous accepted rate, unscaled at block changes",
    "pilot_steps": 50, "series_execution_authorized": False,
    "isolated_level0_wall_s": 7200., "isolated_series_wall_s": 28800.,
    "analysis_stage_wall_s": 3600., "analysis_total_wall_s": 7200.,
    "full_CHNS_study_wall_s": 7200., "required_hashes": HASHES}


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def append_json(path, value):
    with Path(path).open("a") as stream:
        stream.write(json.dumps(value, allow_nan=False)+"\n")


def freeze_policy():
    folder = RESULTS/"policy"; folder.mkdir(parents=True, exist_ok=True)
    record = {"policy": POLICY, "policy_sha256": json_hash(POLICY),
        "design_sha256": hashlib.sha256(Path("STEP3A5_DESIGN.md").read_bytes()).hexdigest()}
    path = folder/"policy.json"
    if path.exists():
        if read_json(path) != record:
            raise ValueError("STOP: frozen STEP3A5 policy/design changed")
    else:
        write_json(path, record)
    return record


def reserve(name):
    global ACTIVE_OUTPUT
    output = RESULTS/name
    if RESULTS.resolve() not in output.resolve().parents:
        raise ValueError("STEP3A5 output only; historical directories read-only")
    policy = freeze_policy()
    plan = read_json(PRIOR/"schedule_optimization/optimized_plan.json")
    if (plan["plan_sha256"] != HASHES["optimized_plan_sha256"] or
            json_hash({k: v for k, v in plan.items() if k != "plan_sha256"}) != plan["plan_sha256"] or
            plan["blocks"][0]["dt"] != TINY_DT):
        raise ValueError("STOP: immutable plan/dt mismatch")
    output.mkdir(parents=True, exist_ok=False); ACTIVE_OUTPUT = output
    write_json(output/"status.json", {"status": "running", "stage": name})
    source = Path(__file__).resolve().parents[1]
    with tarfile.open(output/"multiphase_source.tar.gz", "x:gz") as archive:
        for path in sorted(source.rglob("*.py")):
            archive.add(path, arcname=str(path.relative_to(source)))
    from ..provenance import source_hash
    from petsc4py import PETSc
    import dolfinx, scipy, ufl
    mumps_header = Path("/usr/local/petsc/linux-gnu-real64-32/include/dmumps_c.h")
    version = re.search(r'#define MUMPS_VERSION\s+"([^"]+)"', mumps_header.read_text())
    provenance = {"git_commit": os.environ.get("STEP3_GIT_COMMIT", "unavailable"),
        "docker_image_digest": os.environ.get("STEP3_DOCKER_IMAGE_DIGEST", "unavailable"),
        "multiphase_source_sha256": source_hash(), "old_formulation_version": ABSOLUTE_VERSION,
        "phase_rate_formulation_version": VERSION, **HASHES, **policy,
        "source_archive_sha256": hashlib.sha256((output/"multiphase_source.tar.gz").read_bytes()).hexdigest(),
        "stack": {"dolfinx": dolfinx.__version__, "PETSc": PETSc.Sys.getVersion(),
            "MUMPS": version.group(1) if version else "unavailable", "scipy": scipy.__version__, "ufl": ufl.__version__}}
    write_json(output/"provenance.json", provenance)
    return output


def require(stage, key="qualification_status", value="passed"):
    record = read_json(RESULTS/stage/"status.json")
    if record.get(key) != value:
        raise ValueError("STOP: prerequisite not passed: "+stage)
    return record


def load_input(output):
    from dolfinx import fem
    s, phi_eq, psi, initial, provenance = historical.scientific_input()
    for key in HASHES:
        if key != "optimized_plan_sha256" and provenance[key] != HASHES[key]:
            raise ValueError("STOP: historical input SHA mismatch: "+key)
    c = s.config
    if (c.snes_atol != POLICY["snes_atol"] or c.snes_rtol != POLICY["snes_rtol"] or
            c.snes_max_it != POLICY["snes_max_it"] or c.quadrature_degree != 12 or c.phase_degree != 2):
        raise ValueError("STOP: options/space/quadrature mismatch")
    field = fem.Function(phi_eq.function_space)
    field.x.array[:] = initial; field.x.scatter_forward()
    record = {**read_json(output/"provenance.json"), **provenance,
        "optimized_plan_sha256": HASHES["optimized_plan_sha256"], "SNES_options": snes_options(c)}
    write_json(output/"provenance.json", record)
    return s, phi_eq, psi, field


def assemble(form):
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    vector = fp.assemble_vector(fem.form(form))
    result = vector.array.copy(); vector.destroy()
    return result


class NewtonTrace:
    """Diagnostic-only monitor; raw PETSc bt trace carries line-search trials."""
    def __init__(self, engine, output):
        from petsc4py import PETSc
        self.engine, self.output, self.rows, self.interval = engine, output, [], 1
        self.started = time.perf_counter()
        snes = engine.problem.solver
        # This installs only a viewer/monitor, never changes the bt algorithm.
        opts = PETSc.Options()
        key = snes.getOptionsPrefix()+"snes_linesearch_monitor"
        opts[key] = "ascii:"+str((output/"line_search.txt").resolve())
        snes.setFromOptions()
        del opts[key]
        snes.cancelMonitor(); snes.setMonitor(self)

    def __call__(self, snes, iteration, residual):
        from petsc4py import PETSc
        ksp = snes.getKSP()
        row = {"interval": self.interval, "iteration": iteration, "residual": residual,
            "elapsed_s": time.perf_counter()-self.started,
            "KSP_iterations": ksp.getIterationNumber(), "KSP_reason": ksp.getConvergedReason(),
            "line_search_type": snes.getLineSearch().getType(), "line_search_lambda": None,
            "lambda_note": "pinned petsc4py exposes no getLambda; see raw line_search.txt trial trace",
            "MUMPS_INFOG_1": None, "Jacobian_infinity_norm": None,
            "Jacobian_assembly_cumulative_s": PETSc.Log.Event("SNESJacobianEval").getPerfInfo()["time"]}
        if iteration:
            try:
                row["Jacobian_infinity_norm"] = snes.getJacobian()[0].norm(PETSc.NormType.NORM_INFINITY)
                row["MUMPS_INFOG_1"] = ksp.getPC().getFactorMatrix().getMumpsInfog(1)
            except (PETSc.Error, AttributeError):
                pass
        self.rows.append(row); append_json(self.output/"snes_iterations.jsonl", row)

    def summary(self):
        snes = self.engine.problem.solver
        rows = [r for r in self.rows if r["interval"] == self.interval]
        r0 = rows[0]["residual"] if rows else None
        from ..nonlinear_accuracy import actual_controls
        controls = actual_controls(snes)
        atol, rtol = controls["atol"], controls["rtol"]
        return {"reason": snes.getConvergedReason(), "iterations": snes.getIterationNumber(),
            "initial_residual": r0, "final_residual": snes.getFunctionNorm(),
            "atol_target": atol, "rtol_target": rtol*r0 if r0 is not None else None,
            "effective_target": max(atol, rtol*r0) if r0 is not None else None,
            "actual_controls": controls,
            "snes_type": snes.getType(), "line_search_type": snes.getLineSearch().getType(),
            "tolerances": list(snes.getTolerances()), "accepted": snes.getConvergedReason() > 0}


def physical_sample(s, diagnostics, perf, output):
    from ..diagnostics import interface_resolution
    from ..ch_validation_diagnostics import WallTrace
    with perf.measure("diagnostics"):
        row = diagnostics.measure()
    with perf.measure("interface_resolution"):
        res = interface_resolution(s)
    with perf.measure("wall_topology"):
        crossings = WallTrace(s).crossings()
    row.update(step=s.step_number, dt=s.current_dt or 0.,
        cells_across_transition_certified_min=res.get("cells_across_transition_certified_min", 0.),
        wall_crossings=len(crossings))
    append_json(output/"resolutions.jsonl", res)
    append_json(output/"history.jsonl", row)
    return row


def physical_gates(old, new, audit):
    return {"finite": bool(np.isfinite(list(new.values())).all()),
        "mass": abs(new["mass_error_relative_to_domain"]) <= POLICY["mass_domain"],
        "resolution": new["cells_across_transition_certified_min"] >= POLICY["certified_cells"],
        "topology": new["wall_crossings"] == 2,
        "energy": new["E_total"]-old["E_total"] <= POLICY["energy_growth_J_per_m"],
        "weak_work": abs(audit["weak_work_defect"]) <= POLICY["weak_work_J_per_m"],
        "decomposition": abs(audit["decomposition_roundoff"]) <= POLICY["weak_work_J_per_m"]}


def moderate():
    from ..isolated_ch import IsolatedCH
    from ..phase_rate import IsolatedCHPhaseRate, PhaseMetric
    from ..diagnostics import Diagnostics
    from ..be_work_diagnostics import BEWorkDiagnostics
    from petsc4py import PETSc
    output = reserve("moderate_dt"); perf = Performance(); PETSc.Log.begin()
    values, fields = {}, {}
    try:
        for name, cls in (("absolute", IsolatedCH), ("phase_rate", IsolatedCHPhaseRate)):
            folder = output/name; folder.mkdir()
            with perf.measure("initialization"):
                s, phi_eq, psi, initial = load_input(output)
                engine = cls(s, initial); engine.performance = perf
                trace = NewtonTrace(engine, folder)
                diagnostic, work = Diagnostics(s), BEWorkDiagnostics(s)
                metric = PhaseMetric(s, phi_eq.function_space)
            old = physical_sample(s, diagnostic, perf, folder)
            if name == "phase_rate":
                write_json(folder/"before_solve_scale.json", engine.rate_statistics(MODERATE_DT))
            log = engine.step(MODERATE_DT)
            new = physical_sample(s, diagnostic, perf, folder)
            with perf.measure("BE_work"):
                audit = work.measure(MODERATE_DT, old, new)
            fields[name] = {"phi": s.state.sub(2).collapse().x.array.copy(),
                "mu": s.state.sub(3).collapse().x.array.copy()}
            np.savez_compressed(folder/"physical.npz", **fields[name])
            values[name] = {"old": old, "new": new, "SNES": trace.summary(),
                "log": log, "BE_work": audit, "phi_L2": metric.norm(fields[name]["phi"]),
                "mu_L2": metric.norm(fields[name]["mu"])}
            write_json(folder/"result.json", values[name])
            print(name+" moderate "+json.dumps(values[name]["SNES"]), flush=True)
        a, b = values["absolute"], values["phase_rate"]
        differences = {k+"_L2_relative": metric.norm(fields["phase_rate"][k]-fields["absolute"][k])/
            metric.norm(fields["absolute"][k]) for k in ("phi", "mu")}
        differences.update(energy_absolute=abs(b["new"]["E_total"]-a["new"]["E_total"]),
            mass_absolute=abs(b["new"]["phase_mass"]-a["new"]["phase_mass"]),
            D_CH_relative=abs(b["new"]["CH_dissipation"]/a["new"]["CH_dissipation"]-1))
        gates = {k: value <= (POLICY["moderate_energy_mass_absolute"] if k.endswith("absolute")
                             else POLICY["moderate_phi_mu_D_relative"]) for k, value in differences.items()}
        gates["both_SNES_converged"] = all(v["SNES"]["accepted"] for v in values.values())
        result = {"status": "complete", "qualification_status": "passed" if all(gates.values()) else "failed",
            "dt": MODERATE_DT, "gates": gates, "differences": differences, "runs": values}
        write_json(output/"status.json", result); print(json.dumps({"gates": gates, "differences": differences}), flush=True)
        return result
    finally:
        write_json(output/"timing.json", {**perf.snapshot(), "PETSc": petsc_events()})


def linear_tiny():
    require("moderate_dt")
    from scipy.sparse import bmat
    from ..linearized_ch import PETScSparseFactor
    from petsc4py import PETSc
    output = reserve("tiny_step_linear"); perf = Performance(); PETSc.Log.begin()
    try:
        with perf.measure("initialization"):
            op, q0, record = historical.load_operator()
            op.resolvent_backend, op.performance = "petsc_mumps", perf
        qnew = op.be_step(q0, TINY_DT)
        with perf.measure("assembly"):
            matrix = bmat([[op.M0, op.mobility*op.K], [-TINY_DT*op.H, op.M0]], format="csr")
        with perf.measure("factorization_setup"):
            factor = PETScSparseFactor(matrix)
        with perf.measure("linear_solve"):
            result = factor.solve(np.concatenate((np.zeros(op.size), op.H @ q0)))
        factor.close(); op.clear_be_cache()
        rate, mu = result[:op.size], result[op.size:]
        qrate = q0+TINY_DT*rate  # no mass fix
        errors = {"endpoint_L2_relative": op.norm(qrate-qnew)/op.norm(qnew),
            "increment_L2_relative": op.norm(TINY_DT*rate-(qnew-q0))/op.norm(qnew-q0),
            "mass_domain": abs(float(op.mass @ (qrate-q0)))/op.area}
        gates = {"endpoint": errors["endpoint_L2_relative"] <= POLICY["linear_endpoint_relative"],
            "increment": errors["increment_L2_relative"] <= POLICY["linear_increment_relative"],
            "mass": errors["mass_domain"] <= POLICY["linear_mass_domain"]}
        np.savez_compressed(output/"linear_step.npz", q0=q0, q_new=qnew, q_rate=qrate, phase_rate=rate, mu=mu)
        status = {"status": "complete", "qualification_status": "passed" if all(gates.values()) else "failed",
            "dt": TINY_DT, "errors": errors, "gates": gates, "operator_fingerprint": op.fingerprint,
            "increment_L2": op.norm(qnew-q0), "rate_L2": op.norm(rate),
            "benchmark": "full sparse LINEAR BE only, not nonlinear evidence"}
        write_json(output/"status.json", status); print(json.dumps(status), flush=True)
        return status
    finally:
        write_json(output/"timing.json", {**perf.snapshot(), "PETSc": petsc_events()})


def moderate_postmortem():
    """Read-only diagnostic: reassemble saved states; NO Newton/PDE time step."""
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from scipy.sparse import csr_matrix
    from ..isolated_ch import IsolatedCH
    from ..phase_rate import IsolatedCHPhaseRate, PhaseMetric
    from ..equilibrium import _measures
    from petsc4py import PETSc
    comparison = read_json(RESULTS/"moderate_dt/status.json")
    output = reserve("algebraic_equivalence/moderate_postmortem")
    perf = Performance(); PETSc.Log.begin()
    try:
        with perf.measure("initialization"):
            s, phi_eq, psi, initial = load_input(output)
            absolute = IsolatedCH(s, initial)
            rate = IsolatedCHPhaseRate(s, initial)
            metric = rate.metric
            absolute.dt.value = rate.dt.value = MODERATE_DT
            ap = np.asarray(absolute.state.function_space.sub(0).collapse()[1])
            am = np.asarray(absolute.state.function_space.sub(1).collapse()[1])
            old_phi, _ = ufl.split(absolute.old)
            phi, mu = ufl.split(absolute.state)
            test, _ = ufl.TestFunctions(absolute.state.function_space)
            dx, _ = _measures(s)
        states, arrays = {}, {}
        for name in ("absolute", "phase_rate"):
            with np.load(RESULTS/"moderate_dt"/name/"physical.npz") as data:
                p, m = data["phi"].copy(), data["mu"].copy()
            recovered_rate = (p-initial.x.array)/MODERATE_DT
            absolute.state.x.array[ap], absolute.state.x.array[am] = p, m
            absolute.state.x.scatter_forward()
            rate.state.x.array[rate.rate_map], rate.state.x.array[rate.mu_map] = recovered_rate, m
            rate.state.x.scatter_forward()
            with perf.measure("assembly"):
                phase = assemble(absolute.phase_form)[ap]
                chemical = assemble(absolute.chemical_form)[am]
                retained_expression = assemble(rate.phase_form)[rate.rate_map]
                rate_chemical = assemble(rate.chemical_form)[rate.mu_map]
                coefficient = fem.Function(phi_eq.function_space)
                coefficient.x.array[:] = recovered_rate; coefficient.x.scatter_forward()
                flux = assemble(s.config.mobility*ufl.inner(ufl.grad(mu), ufl.grad(test))*dx)[ap]
                coefficient_ufl = assemble(coefficient*ufl.TestFunction(phi_eq.function_space)*dx)+flux
                sparse = metric.M @ recovered_rate+s.config.mobility*(metric.K @ m)
                time_literal = assemble((phi-old_phi)/absolute.dt*test*dx)[ap]
            residual_vectors = {"absolute_literal_phase": phase, "physical_chemical": chemical,
                "coefficient_difference_UFL_phase": coefficient_ufl,
                "recovered_rate_UFL_phase": retained_expression, "sparse_coefficient_phase": sparse,
                "recovered_rate_chemical": rate_chemical}
            states[name] = {"algebraic_norms": {k: float(np.linalg.norm(v)) for k, v in residual_vectors.items()},
                "Riesz_L2_norms": {k: metric.riesz_norm(v) for k, v in residual_vectors.items()},
                "chemical_form_difference_norm": float(np.linalg.norm(rate_chemical-chemical)),
                "literal_time_vs_coefficient_norm": float(np.linalg.norm(time_literal-metric.M @ recovered_rate)),
                "phase_time_block_norm": float(np.linalg.norm(time_literal)),
                "phase_flux_block_norm": float(np.linalg.norm(flux)),
                "rate_recovered_L2": metric.norm(recovered_rate),
                "rate_recovered_max": float(max(abs(recovered_rate))),
                "dt_times_rate_recovered_max": float(MODERATE_DT*max(abs(recovered_rate))),
                "max_ULP_over_dt": float(max(abs(np.spacing(p)))/MODERATE_DT),
                "rate_recovery_note": "diagnostic coefficient difference only; original moderate rate coefficients were not archived",
                "reconstruction_bitwise_equal": bool(np.array_equal(initial.x.array+MODERATE_DT*recovered_rate, p))}
            arrays[name] = {"p": p, "m": m, "phase": phase, "chemical": chemical}
        # A Jacobian MATVEC explains the measured displacement; it is NOT a
        # Newton update, not a solve and is never published into physical state.
        absolute.state.x.array[ap] = arrays["absolute"]["p"]
        absolute.state.x.array[am] = arrays["absolute"]["m"]
        absolute.state.x.scatter_forward()
        Jform = ufl.derivative(absolute.F, absolute.state, ufl.TrialFunction(absolute.state.function_space))
        with perf.measure("assembly"):
            matrix = fp.assemble_matrix(fem.form(Jform)); matrix.assemble()
            ptr, col, val = matrix.getValuesCSR()
            J = csr_matrix((val.copy(), col.copy(), ptr.copy()), shape=matrix.getSize()); matrix.destroy()
        difference = np.zeros(J.shape[0])
        difference[ap] = arrays["phase_rate"]["p"]-arrays["absolute"]["p"]
        difference[am] = arrays["phase_rate"]["m"]-arrays["absolute"]["m"]
        action = J @ difference
        closure = arrays["absolute"]["chemical"]+action[am]-arrays["phase_rate"]["chemical"]
        dmu = difference[am]
        mean = float(metric.mass @ dmu/metric.area)
        result = {"status": "complete", "qualification_status": "diagnostic_only",
            "PDE_steps": 0, "nonlinear_solves": 0, "no_tolerance_or_policy_changes": True,
            "comparison_qualification": comparison["qualification_status"], "states": states,
            "chemical_Jacobian_displacement_closure_relative": float(np.linalg.norm(closure)/
                np.linalg.norm(arrays["absolute"]["chemical"])),
            "mu_difference_L2": metric.norm(dmu), "mu_difference_mean": mean,
            "mu_difference_zero_mean_L2": metric.norm(dmu-mean),
            "source_comparison_sha256": hashlib.sha256((RESULTS/"moderate_dt/status.json").read_bytes()).hexdigest(),
            "interpretation": "Both original solvers stopped at their unchanged SNES targets; this audit does not override the failed physical mu/D comparison gate."}
        write_json(output/"status.json", result); print(json.dumps(result, indent=2), flush=True)
        return result
    finally:
        write_json(output/"timing.json", {**perf.snapshot(), "PETSc": petsc_events()})


STAGES = {"freeze_policy": freeze_policy, "moderate": moderate, "linear_tiny": linear_tiny,
    "moderate_postmortem": moderate_postmortem}
