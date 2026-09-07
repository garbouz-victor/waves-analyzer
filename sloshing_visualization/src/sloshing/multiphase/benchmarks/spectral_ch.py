"""STEP 3A.3 staged analysis. Historical inputs are always read-only."""
import csv
import json
import hashlib
import os
from pathlib import Path
import tarfile
import time
import numpy as np
from ..linearized_ch import array_hash, json_hash

RESULTS = Path("validation_results/step3a3")
EQUILIBRIUM = Path("validation_results/step3a2/equilibrium60/prepared_symmetric")
OLD_BUMP = Path("validation_results/step3a2/equilibrium_perturbation/be_dt_1e-5")
ACTIVE_OUTPUT = None


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def reserve(name):
    global ACTIVE_OUTPUT
    output = RESULTS/name
    if RESULTS.resolve() not in output.resolve().parents:
        raise ValueError("Historical/external output paths are forbidden; require a STEP3A3 subdirectory")
    output.mkdir(parents=True, exist_ok=False)
    ACTIVE_OUTPUT = output
    write_json(output/"status.json", {"status": "running", "stage": name})
    source = Path(__file__).resolve().parents[1]
    with tarfile.open(output/"multiphase_source.tar.gz", "x:gz") as archive:
        for path in sorted(source.rglob("*.py")):
            archive.add(path, arcname=str(path.relative_to(source)))
    from ..provenance import source_hash
    inputs = RESULTS/"linear_operator/status.json"
    provenance = read_json(inputs).get("provenance", {}) if inputs.exists() else {}
    provenance.update(git_commit=os.environ.get("STEP3_GIT_COMMIT", "unavailable"),
        docker_image_digest=os.environ.get("STEP3_DOCKER_IMAGE_DIGEST", "unavailable"),
        multiphase_source_sha256=source_hash(),
        source_archive_sha256=hashlib.sha256((output/"multiphase_source.tar.gz").read_bytes()).hexdigest(),
        design_sha256=hashlib.sha256(Path("STEP3A3_DESIGN.md").read_bytes()).hexdigest())
    write_json(output/"provenance.json", provenance)
    return output


def scientific_input():
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..equilibrium import load_prepared, zero_mass_perturbations
    from ..provenance import run_provenance
    meta = read_json(EQUILIBRIUM/"prepared.json")
    s = CHNSSolver(ModelConfig(**meta["config"]))
    eq = load_prepared(EQUILIBRIUM, s, meta["target_mass"])
    s.initialize_prepared_equilibrium(eq)
    phi = s.state.sub(2).collapse()
    old = read_json(OLD_BUMP/"perturbation.json")
    if old["delta"] != 1e-3:
        raise ValueError("Historical amplitude changed")
    frozen = RESULTS/"linear_operator/inputs.npz"
    if frozen.exists():
        with np.load(frozen) as data:
            psi = data["psi"].copy()
            initial = data["initial_phi"].copy()
            if not np.array_equal(data["phi_eq"], phi.x.array):
                raise ValueError("Frozen equilibrium bytes mismatch")
        recovery = "load frozen, historical-hash-verified coefficient arrays"
    else:
        field, _ = zero_mass_perturbations(s, phi)["smooth_bulk"]
        psi = field.x.array.copy()
        initial = phi.x.array+1e-3*psi
        recovery = "replay unchanged historical initializer; require BOTH archived coefficient hashes"
    if (array_hash(psi) != old["psi_coefficient_sha256"] or
            array_hash(initial) != old["initial_phi_coefficient_sha256"] or
            not np.array_equal(initial, phi.x.array+1e-3*psi)):
        raise ValueError("STOP: perturbation/initial phi not coefficientwise identical to STEP 3A.2")
    provenance = {**run_provenance(s), "equilibrium_fingerprint": meta["fingerprint"],
        "perturbation_coefficient_sha256": array_hash(psi), "initial_phi_sha256": array_hash(initial),
        "quadrature_degree": s.config.quadrature_degree, "mobility": s.config.mobility,
        "target_mass": meta["target_mass"], "input_recovery": recovery}
    return s, phi, psi, initial, provenance


def historical_fast_balance():
    rows = []
    times_per_step = []
    root = OLD_BUMP.parent
    for name in ("be_dt_1e-5", "be_dt_5e-6", "be_dt_2_5e-6"):
        summary = read_json(root/name/"summary.json")
        times_per_step.append(summary["runtime_s"]/summary["accepted_steps"])
        with (root/name/"history.csv").open() as stream:
            old = list(csv.DictReader(stream))
        for row in old[:4]:
            values = {k: float(row[k]) for k in ("time", "CH_dissipation", "viscous_dissipation", "slip_dissipation", "E_kin")}
            values.update(run=name,
                hydrodynamic_power_over_CH=(values["viscous_dissipation"]+values["slip_dissipation"])/values["CH_dissipation"],
                kinetic_over_initial_excess=values["E_kin"]/5.542888467657825e-8)
            rows.append(values)
    return {"measurements": rows, "max_hydrodynamic_power_over_CH": max(r["hydrodynamic_power_over_CH"] for r in rows),
        "max_kinetic_over_initial_excess": max(r["kinetic_over_initial_excess"] for r in rows),
        "historical_seconds_per_step": times_per_step,
        "cost_seconds_per_step_conservative": max(times_per_step),
        "scope": "CH-dominated stiff subsystem planner; final qualification remains FULL nonlinear CHNS"}


def operator_stage():
    from scipy.sparse import save_npz
    from ..linearized_ch import assemble_linearized_ch, NonlinearCHRate, ALGORITHM_VERSION
    output = reserve("linear_operator")
    started = time.perf_counter()
    try:
        s, phi, psi, initial, provenance = scientific_input()
        np.savez_compressed(output/"inputs.npz", phi_eq=phi.x.array, psi=psi, initial_phi=initial)
        provenance["spectrum_algorithm_version"] = ALGORITHM_VERSION
        op = assemble_linearized_ch(s, phi, provenance)
        checks = op.matrix_checks()
        nonlinear = NonlinearCHRate(s, phi.function_space, op)
        finite = nonlinear.finite_difference_checks(phi.x.array, psi)
        # Additional independent FE directions, not only the seed.
        rng = np.random.default_rng(3302)
        random_finite = []
        for _ in range(3):
            q = op.project(rng.normal(size=op.size)); q /= max(abs(q))
            random_finite.append(nonlinear.finite_difference_checks(phi.x.array, q))
        balance = historical_fast_balance()
        q0 = op.project(psi)
        timescale = op.energy(q0)/op.power(q0)
        gates = {"mass_conservation": checks["max_relative_mass_defect"] <= 1e-11,
            "roundoff_projection": checks["max_projection_relative_defect"] <= 1e-11,
            "matrix_symmetry": max(checks["frobenius_symmetry_defects"].values()) <= 1e-11 and
                               max(checks["bilinear_symmetry_defects"].values()) <= 1e-11,
            "Hessian_linearization": all(r[-1]["hessian_relative_M_dual_error"] <= 1e-6 for r in [finite]+random_finite),
            "nonlinear_rate_linearization": all(r[-1]["central_rate_relative_L2_error"] <= 1e-6 for r in [finite]+random_finite),
            "linear_energy_identity": checks["max_relative_energy_identity_defect"] <= 1e-10,
            "CH_dominated": balance["max_hydrodynamic_power_over_CH"] <= 1e-4 and balance["max_kinetic_over_initial_excess"] <= 1e-4,
            "initial_timescale_agreement": abs(timescale/4.862048987156724e-7-1) <= 1e-2}
        for name in ("M0", "K", "H"):
            save_npz(output/(name+".npz"), getattr(op, name))
        np.savez_compressed(output/"mass.npz", mass=op.mass)
        write_json(output/"operator_identity.json", {**op.identity, "fingerprint": op.fingerprint})
        from petsc4py import PETSc
        import scipy, dolfinx
        import importlib.util
        stack = {"dolfinx": dolfinx.__version__, "petsc": PETSc.Sys.getVersion(), "scipy": scipy.__version__,
                 "slepc4py_available": importlib.util.find_spec("slepc4py") is not None}
        if stack["slepc4py_available"]:
            from slepc4py import SLEPc
            stack["slepc"] = SLEPc.Sys.getVersion()
        result = {"status": "complete", "qualification_status": "passed" if all(gates.values()) else "failed",
            "gates": gates, "provenance": provenance, "operator_fingerprint": op.fingerprint,
            "stack": stack, "phase_dofs": op.size, "matrix_checks": checks,
            "smooth_bump_finite_differences": finite, "random_finite_differences": random_finite,
            "CH_dominated_justification": balance, "initial_linear_E2_at_delta": 1e-6*op.energy(q0),
            "initial_linear_D2_at_delta": 1e-6*op.power(q0), "tau_E_linear": timescale,
            "historical_tau_E": 4.862048987156724e-7,
            "old_dt_over_tau_E": [dt/4.862048987156724e-7 for dt in (1e-5, 5e-6, 2.5e-6)],
            "raw_psi_mass": float(op.mass @ psi), "projected_psi_mass": float(op.mass @ q0),
            "runtime_s": time.perf_counter()-started}
        write_json(output/"status.json", result)
        print(json.dumps(result, indent=2), flush=True)
        return result
    except BaseException as error:
        write_json(output/"failure.json", {"status": "failed", "error": str(error), "runtime_s": time.perf_counter()-started})
        raise


def load_operator():
    from scipy.sparse import load_npz
    from ..linearized_ch import LinearizedCH
    folder = RESULTS/"linear_operator"
    result = read_json(folder/"status.json")
    if result.get("qualification_status") != "passed":
        raise ValueError("STOP: linear operator prerequisite failed")
    with np.load(folder/"mass.npz") as data:
        mass = data["mass"].copy()
    op = LinearizedCH(*(load_npz(folder/(name+".npz")) for name in ("M0", "K", "H")),
        mobility=result["provenance"]["mobility"], mass_vector=mass, provenance=result["provenance"])
    if op.fingerprint != result["operator_fingerprint"]:
        raise ValueError("Saved operator fingerprint mismatch")
    with np.load(folder/"inputs.npz") as data:
        psi = data["psi"].copy()
    if array_hash(psi) != result["provenance"]["perturbation_coefficient_sha256"]:
        raise ValueError("Saved perturbation hash mismatch")
    return op, op.project(psi), result


def spectrum_stage():
    from ..ch_krylov import MassArnoldi, compare_references, ALGORITHM_VERSION
    from ..provenance import source_hash
    output = reserve("spectrum")
    krylov_output = RESULTS/"krylov"
    krylov_output.mkdir(exist_ok=False)
    started = time.perf_counter()
    op, q0, operator_record = load_operator()
    provenance = {**operator_record["provenance"], "multiphase_source_sha256": source_hash(),
                  "spectrum_algorithm_version": ALGORITHM_VERSION, "operator_fingerprint": op.fingerprint}
    levels = (32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024)
    basis_bytes = 2*op.size*(levels[-1]+1)*8
    if basis_bytes > 2*1024**3:
        raise RuntimeError("STOP: Arnoldi basis resource guard exceeded")
    arnoldi = MassArnoldi(op, q0, levels[-1])
    previous, previous_ritz, records = None, None, []
    qualified = False
    reason = "Krylov self-convergence not reached at finite resource guard"
    times = None
    for m in levels:
        if time.perf_counter()-started > 3600:
            reason = "Krylov one-hour investigation guard exceeded"
            break
        def monitor(j):
            if j % 16 == 0:
                print(f"Arnoldi dimension {j}; elapsed {time.perf_counter()-started:.2f}s", flush=True)
        arnoldi.extend(m, monitor)
        ritz = arnoldi.ritz(previous=previous_ritz, dissipative_significance=True)
        if times is None:
            rho_estimate = max(abs(r["real"]) for r in ritz)
            times = np.concatenate(([0.], np.geomspace(1e-3/rho_estimate, 1e-4, 81)))
        reference = arnoldi.reference(times)
        checks = arnoldi.diagnostics()
        comparison = compare_references(op, previous, reference) if previous is not None else None
        unstable = [r for r in ritz if r["significant"] and r["converged"] and
                    r["real"] > 1e-8*max(np.hypot(r["real"], r["imag"]), 1.)]
        complex_modes = [r for r in ritz if r["significant"] and r["converged"] and
                         abs(r["imag"])/max(abs(r["real"]), 1.) > 1e-8]
        self_pass = comparison is not None and max(comparison[k] for k in (
            "max_relative_L2_difference", "max_relative_energy_difference", "max_relative_power_difference")) <= 1e-6
        relevant = [r for r in ritz if r["converged"] and r["significant"] and r["real"] < 0]
        gates = {"orthogonality": checks["orthogonality_defect"] <= 1e-11,
                 "no_significant_positive_eigenvalues": not unstable,
                 "no_significant_complex_eigenvalues": not complex_modes,
                 "reference_self_convergence": self_pass, "converged_relevant_modes_present": bool(relevant),
                 "minimum_three_levels": len(records) >= 2}
        record = {"dimension": arnoldi.dimension, "diagnostics": checks, "gates": gates,
            "reference_comparison": comparison, "ritz": ritz,
            "max_projected_energy_identity_relative_defect": float(max(reference["projected_energy_identity_relative_defect"])),
            "runtime_s": time.perf_counter()-started, "provenance": provenance}
        records.append(record)
        write_json(output/f"level_m{m}.json", record)
        print(json.dumps({k: record[k] for k in ("dimension", "gates", "reference_comparison", "runtime_s")}), flush=True)
        previous, previous_ritz = reference, ritz
        if unstable or complex_modes or not gates["orthogonality"]:
            reason = "STOP: significant spectral sign/imaginary or orthogonality failure"
            break
        if all(gates.values()):
            qualified = True
            reason = "All reference/spectrum prerequisites passed"
            break
    m = arnoldi.dimension
    np.savez_compressed(krylov_output/"basis.npz", V=arnoldi.V[:, :m+1],
                        hessenberg=arnoldi.hessenberg[:m+1, :m], seed=q0, beta=arnoldi.beta)
    np.savez_compressed(krylov_output/"krylov_reference.npz", **previous)
    fingerprint = arnoldi.fingerprint()
    relevant = [r for r in previous_ritz if r["converged"] and r["significant"] and r["real"] < 0]
    rho = max(abs(r["real"]) for r in relevant) if relevant else None
    with (output/"ritz.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=previous_ritz[0].keys())
        writer.writeheader(); writer.writerows(previous_ritz)
    result = {"status": "complete", "qualification_status": "passed" if qualified else "failed",
        "reason": reason, "provenance": provenance, "operator_fingerprint": op.fingerprint,
        "Krylov_sha256": fingerprint, "dimension": m, "levels": [r["dimension"] for r in records],
        "last_reference_comparison": records[-1]["reference_comparison"],
        "rho_relevant": rho, "tau_fast_relevant": 1/rho if rho else None,
        "last_gates": records[-1]["gates"], "ritz": previous_ritz, "basis_storage_bytes": basis_bytes,
        "reference_sha256": array_hash(previous["states"]), "runtime_s": time.perf_counter()-started,
        "diagnostic_only": True, "nonlinear_qualification": "not run"}
    write_json(output/"spectrum.json", result)
    write_json(output/"status.json", result)
    write_json(krylov_output/"status.json", result)
    print(json.dumps({k: result[k] for k in ("qualification_status", "reason", "dimension", "rho_relevant", "tau_fast_relevant", "runtime_s")}), flush=True)
    return result


def restarted_reference_stage():
    from ..ch_krylov import RestartedExponential, compare_references
    from ..provenance import source_hash
    output = reserve("krylov/restarted")
    op, q0, operator_record = load_operator()
    spectrum = read_json(RESULTS/"spectrum/spectrum.json")
    for key in ("orthogonality", "no_significant_positive_eigenvalues", "no_significant_complex_eigenvalues"):
        if not spectrum["last_gates"][key]:
            raise ValueError("STOP: unresolved spectral sign/orthogonality failure")
    with np.load(RESULTS/"krylov/krylov_reference.npz") as data:
        times = data["times"].copy()
    provenance = {**operator_record["provenance"], "multiphase_source_sha256": source_hash(),
        "spectrum_algorithm_version": "SLEPc-MFN-restarted-exponential-M-certified-v1",
        "operator_fingerprint": op.fingerprint,
        "internal_inner_product": "Euclidean SLEPc MFN; external self-convergence in FEM M norm"}
    started = time.perf_counter()
    previous, records = None, []
    passed = False
    for ncv in (32, 48, 64, 96, 128):
        if time.perf_counter()-started > 3600:
            break
        engine = RestartedExponential(op, dimension=ncv, tolerance=1e-11)
        def monitor(row):
            if row["operator_applications"] > 100 or row["time"] == 0:
                print(f"Restarted ncv={ncv} t={row['time']:.8g} applications={row['operator_applications']} "
                      f"iterations={row['iterations']} elapsed={time.perf_counter()-started:.2f}s", flush=True)
        try:
            reference, history = engine.reference(q0, times, monitor)
        finally:
            engine.close()
        comparison = compare_references(op, previous, reference) if previous is not None else None
        converged = comparison is not None and max(comparison[k] for k in (
            "max_relative_L2_difference", "max_relative_energy_difference", "max_relative_power_difference")) <= 1e-6
        record = {"ncv": ncv, "internal_tolerance": 1e-11, "history": history,
            "reference_comparison": comparison, "self_converged": bool(converged),
            "runtime_s": time.perf_counter()-started, "provenance": provenance,
            "maximum_state_mass": float(max(abs(op.mass @ reference["states"].T)))}
        write_json(output/f"level_ncv{ncv}.json", record)
        np.savez_compressed(output/f"reference_ncv{ncv}.npz", **reference)
        records.append(record)
        previous = reference
        print(json.dumps({k: record[k] for k in ("ncv", "self_converged", "reference_comparison", "runtime_s")}), flush=True)
        if converged and len(records) >= 3:
            passed = True
            break
    result = {"status": "complete", "qualification_status": "passed" if passed else "failed",
        "reason": "Restarted exponential self-converged" if passed else "STOP: restarted reference did not pass finite resource guard",
        "levels": [r["ncv"] for r in records], "last_reference_comparison": records[-1]["reference_comparison"],
        "provenance": provenance, "operator_fingerprint": op.fingerprint,
        "rho_relevant": spectrum["rho_relevant"], "tau_fast_relevant": spectrum["tau_fast_relevant"],
        "reference_sha256": array_hash(previous["states"]),
        "unrestarted_reference_status": "failed, retained separately; never substituted as a qualified reference",
        "runtime_s": time.perf_counter()-started, "diagnostic_only": True, "nonlinear_qualification": "not run"}
    result["Krylov_sha256"] = json_hash(result)
    write_json(output/"status.json", result)
    return result


def rational_reference_stage(name="krylov/rational_balanced"):
    from ..ch_krylov import RationalMassKrylov, compare_references
    from ..provenance import source_hash
    output = reserve(name)
    op, q0, operator_record = load_operator()
    op.resolvent_backend = "petsc_mumps"
    spectrum = read_json(RESULTS/"spectrum/spectrum.json")
    for key in ("orthogonality", "no_significant_positive_eigenvalues", "no_significant_complex_eigenvalues"):
        if not spectrum["last_gates"][key]:
            raise ValueError("STOP: unresolved primary spectral sign/orthogonality failure")
    rho = spectrum["rho_relevant"]
    if not rho:
        raise ValueError("No converged relevant primary Ritz spectrum")
    shifts = np.geomspace(1/rho, 1e-4, 8)
    with np.load(RESULTS/"krylov/krylov_reference.npz") as data:
        times = data["times"].copy()
    provenance = {**operator_record["provenance"], "multiphase_source_sha256": source_hash(),
        "spectrum_algorithm_version": "fixed-positive-pole-rational-M-Krylov-v1", "operator_fingerprint": op.fingerprint,
        "resolvent_auxiliary_scaling": "r_tilde=sqrt(resolvent_time)*r; same eliminated operator",
        "resolvent_backend": "pinned PETSc MUMPS; at most three simultaneous cached factors"}
    arnoldi = RationalMassKrylov(op, q0, shifts, max_dimension=128)
    previous, previous_ritz, records = None, None, []
    started = time.perf_counter()
    passed = False
    for m in (32, 48, 64, 96, 128):
        if time.perf_counter()-started > 3600:
            break
        def monitor(j):
            if j % 8 == 0:
                print(f"Rational M-Krylov dimension {j}; elapsed {time.perf_counter()-started:.2f}s", flush=True)
        arnoldi.extend(m, monitor)
        ritz = arnoldi.ritz(previous=previous_ritz)
        reference = arnoldi.reference(times)
        comparison = compare_references(op, previous, reference) if previous is not None else None
        checks = arnoldi.diagnostics()
        bad_sign = [r for r in ritz if r["significant"] and r["converged"] and r["real"] > 1e-8*max(abs(r["real"]), 1.)]
        bad_imag = [r for r in ritz if r["significant"] and r["converged"] and abs(r["imag"])/max(abs(r["real"]), 1.) > 1e-8]
        self_pass = comparison is not None and max(comparison[k] for k in (
            "max_relative_L2_difference", "max_relative_energy_difference", "max_relative_power_difference")) <= 1e-6
        gates = {"reference_self_convergence": self_pass, "orthogonality": checks["orthogonality_defect"] <= 1e-11,
                 "no_significant_positive_eigenvalues": not bad_sign,
                 "no_significant_complex_eigenvalues": not bad_imag, "minimum_three_levels": len(records) >= 2}
        record = {"dimension": arnoldi.dimension, "diagnostics": checks, "gates": gates,
            "reference_comparison": comparison, "ritz": ritz, "provenance": provenance,
            "max_projected_energy_identity_relative_defect": float(max(reference["projected_energy_identity_relative_defect"])),
            "runtime_s": time.perf_counter()-started}
        write_json(output/f"level_m{m}.json", record)
        records.append(record)
        print(json.dumps({"dimension": m, "gates": gates,
            "errors": {k: v for k, v in (comparison or {}).items() if k != "relative_L2_by_time"},
            "runtime_s": record["runtime_s"]}), flush=True)
        previous, previous_ritz = reference, ritz
        if bad_sign or bad_imag or not gates["orthogonality"]:
            break
        if all(gates.values()):
            passed = True
            break
    op.clear_be_cache()
    m = arnoldi.dimension
    np.savez_compressed(output/"basis.npz", V=arnoldi.V[:, :m], projected_L=arnoldi.projected(),
                        seed=q0, beta=arnoldi.beta, shifts=shifts)
    np.savez_compressed(output/"krylov_reference.npz", **previous)
    result = {"status": "complete", "qualification_status": "passed" if passed else "failed",
        "reason": "Rational M-Krylov reference passed" if passed else "STOP: rational reference prerequisite failed",
        "levels": [r["dimension"] for r in records], "dimension": m,
        "last_reference_comparison": records[-1]["reference_comparison"], "last_gates": records[-1]["gates"],
        "provenance": provenance, "operator_fingerprint": op.fingerprint,
        "rho_relevant": rho, "tau_fast_relevant": 1/rho, "resolvent_times": shifts.tolist(),
        "reference_sha256": array_hash(previous["states"]), "Krylov_sha256": arnoldi.fingerprint(),
        "ritz": previous_ritz, "runtime_s": time.perf_counter()-started, "diagnostic_only": True}
    write_json(output/"status.json", result)
    return result


def qualified_reference():
    folder = RESULTS/"krylov/rational_balanced"
    result = read_json(folder/"status.json")
    op, seed, operator_record = load_operator()
    if (result.get("qualification_status") != "passed" or
            result["operator_fingerprint"] != op.fingerprint):
        raise ValueError("STOP: no matching qualified Krylov reference")
    with np.load(folder/"basis.npz") as data:
        basis = {key: data[key].copy() for key in data.files}
    with np.load(folder/"krylov_reference.npz") as data:
        reference = {key: data[key].copy() for key in data.files}
    last = read_json(folder/f"level_m{result['dimension']}.json")
    if (array_hash(basis["V"]) != last["diagnostics"]["basis_sha256"] or
            array_hash(basis["projected_L"]) != last["diagnostics"]["projected_L_sha256"]):
        raise ValueError("Krylov basis/projected-operator bytes changed")
    if last["max_projected_energy_identity_relative_defect"] > 1e-10:
        raise ValueError("STOP: reference projected energy identity not qualified")
    if array_hash(reference["states"]) != result["reference_sha256"]:
        raise ValueError("Krylov reference bytes changed")
    if not np.array_equal(basis["seed"], seed):
        raise ValueError("Krylov seed changed")
    return op, seed, operator_record, result, basis, reference


def low_mode_stage():
    from dolfinx import fem
    from ..linearized_ch import shifted_modes
    from ..isolated_ch import IsolatedCH
    from ..diagnostics import Diagnostics, interface_resolution
    from ..be_work_diagnostics import BEWorkDiagnostics
    from .contact_angle import observe
    op, _, _, reference_record, _, _ = qualified_reference()
    output = reserve("low_mode")
    started = time.perf_counter()
    rows, modes = shifted_modes(op, -1/3e-4)
    write_json(output/"eigenmodes.json", rows)
    eligible = [i for i, r in enumerate(rows) if r["real"] < 0 and
        abs(r["imag"])/abs(r["real"]) <= 1e-8 and r["relative_full_operator_residual"] <= 1e-8
        and 1e-4 <= r["tau"] <= 1e-3]
    if not eligible:
        raise RuntimeError("STOP: no fully residual-certified low mode in the predeclared time range")
    j = min(eligible, key=lambda i: abs(np.log(rows[i]["tau"]/3e-4)))
    mode = op.project(modes[:, j].real); mode /= op.norm(mode)
    if mode[np.argmax(abs(mode))] < 0:
        mode = -mode
    rate, tau = rows[j]["real"], rows[j]["tau"]
    mode_record = {**rows[j], "coefficient_sha256": array_hash(mode), "M_norm": op.norm(mode),
        "coefficient_max": float(max(abs(mode))), "operator_fingerprint": op.fingerprint,
        "Krylov_sha256": reference_record["Krylov_sha256"], "target_tau": 3e-4,
        "initial_max_delta_phi": 1e-4, "amplitude_checks": [5e-5, 2.5e-5],
        "temporal_steps": [20, 40, 80], "amplitude_check_steps": 20,
        "benchmark": "isolated CH temporal benchmark; u=0; not full CHNS qualification"}
    write_json(output/"selected_mode.json", mode_record)
    np.savez_compressed(output/"mode.npz", coefficients=mode)
    s, phi_eq, _, _, provenance = scientific_input()
    experiments, temporal, amplitude_checks = [], [], []
    for steps, max_delta in ((20, 1e-4), (40, 1e-4), (80, 1e-4), (20, 5e-5), (20, 2.5e-5)):
        folder = output/f"be_n{steps}_maxdelta_{max_delta:.8g}"
        folder.mkdir(exist_ok=False)
        amplitude = max_delta/max(abs(mode))
        initial = fem.Function(phi_eq.function_space)
        initial.x.array[:] = phi_eq.x.array+amplitude*mode
        initial.x.scatter_forward()
        engine = IsolatedCH(s, initial)
        diagnostics, work = Diagnostics(s), BEWorkDiagnostics(s)
        dt = tau/steps
        histories, audits, snes_history = [], [], []
        write_json(folder/"provenance.json", {**provenance, "selected_mode": mode_record,
            "initial_phi_sha256": array_hash(initial.x.array), "delta_mode": float(amplitude),
            "dt": dt, "steps": steps, "t_end": tau})
        first_geometry = observe(s)
        if len(first_geometry["crossings"]) != 2:
            raise RuntimeError("STOP: low-mode initial contact topology changed")
        def sample():
            row = diagnostics.measure()
            resolution = interface_resolution(s)
            current = s.state.sub(2).collapse().x.array
            q = (current-phi_eq.x.array)/amplitude
            row["modal_amplitude"] = float(op.inner(mode, q).real)
            row["linear_BE_amplitude"] = float((1-rate*dt)**(-s.step_number))
            row["exponential_amplitude"] = float(np.exp(rate*s.time))
            row["cells_across_transition_certified_min"] = resolution.get("cells_across_transition_certified_min", 0.)
            row["mass_error_to_target_domain"] = float(abs(row["phase_mass"]-provenance["target_mass"])/op.area)
            with (folder/"resolutions.jsonl").open("a") as stream:
                stream.write(json.dumps(resolution)+"\n")
            histories.append(row)
            if row["mass_error_to_target_domain"] > 1e-10 or not resolution["qualified_resolution"]:
                raise RuntimeError("STOP: low-mode mass/spatial gate failure")
            return row
        first = previous = sample()
        with (folder/"history.csv").open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=first.keys()); writer.writeheader(); writer.writerow(first)
            for n in range(steps):
                log = engine.step(dt, (n+1)*dt)
                row = sample()
                audit = work.measure(dt, previous, row); audit.update(interval=n+1, time=s.time)
                audits.append(audit); snes_history.append(log)
                writer.writerow(row); stream.flush()
                with (folder/"be_work.jsonl").open("a") as audit_stream:
                    audit_stream.write(json.dumps(audit)+"\n")
                with (folder/"snes.jsonl").open("a") as snes_stream:
                    snes_stream.write(json.dumps(log)+"\n")
                if row["E_total"]-previous["E_total"] > 1e-9:
                    raise RuntimeError("STOP: low-mode energy growth")
                previous = row
                if n == 0 or (n+1) % 10 == 0:
                    print(f"Low mode n={steps} amplitude={max_delta:g}: {n+1}/{steps}, "
                          f"modal={row['modal_amplitude']:.10g} D_CH={row['CH_dissipation']:.6g}", flush=True)
        end = histories[-1]
        change = end["E_total"]-first["E_total"]
        measured = end["modal_amplitude"]
        expected_be = end["linear_BE_amplitude"]
        final_geometry = observe(s)
        experiment = {"steps": steps, "dt": dt, "t_end": tau, "max_delta_phi": max_delta,
            "delta_mode": float(amplitude), "measured_modal_amplitude": measured,
            "linear_BE_amplitude": expected_be, "exponential_amplitude": float(np.exp(-1)),
            "error_to_exponential": float(abs(measured-np.exp(-1))),
            "relative_error_to_BE": float(abs(measured/expected_be-1)),
            "DeltaE": change, "integrated_D_CH": end["cumulative_CH_dissipation"],
            "final_budget_defect": end["energy_budget_defect"],
            "final_budget_over_actual_change": abs(end["energy_budget_defect"])/max(abs(change), 1e-300),
            "max_mass_domain_error": max(r["mass_error_to_target_domain"] for r in histories),
            "minimum_certified_cells": min(r["cells_across_transition_certified_min"] for r in histories),
            "max_energy_growth": max(r["delta_E_interval"] for r in histories),
            "max_weak_work_defect": max(abs(r["weak_work_defect"]) for r in audits),
            "first_BE_work": audits[0], "runtime_s": sum(r["runtime_s"] for r in snes_history),
            "initial_phi_sha256": array_hash(initial.x.array),
            "geometry": {"initial": first_geometry, "final": final_geometry}}
        write_json(folder/"summary.json", experiment)
        np.savez_compressed(folder/"endpoint.npz", phi=s.state.sub(2).collapse().x.array)
        experiments.append(experiment)
        (temporal if max_delta == 1e-4 else amplitude_checks).append(experiment)
    errors = [r["error_to_exponential"] for r in temporal]
    powers = [r["integrated_D_CH"] for r in temporal]
    power_diff = abs(np.diff(powers))
    amp_errors = [temporal[0]["relative_error_to_BE"]]+[r["relative_error_to_BE"] for r in amplitude_checks]
    orders = np.log2(np.asarray(errors[:-1])/errors[1:])
    gates = {"mass": all(r["max_mass_domain_error"] <= 1e-10 for r in experiments),
        "resolution": all(r["minimum_certified_cells"] >= 8 for r in experiments),
        "correct_decay": all(0 < r["measured_modal_amplitude"] < 1 for r in experiments),
        "first_order_decay": bool(min(orders) >= .8),
        "linear_BE_small_amplitude": amp_errors[-1] <= 1e-3,
        "amplitude_limit": all(b <= max(a, 1e-8) for a, b in zip(amp_errors, amp_errors[1:])),
        "energy_decreases": all(r["max_energy_growth"] <= 1e-9 and r["DeltaE"] < 0 for r in experiments),
        "integrated_D_converges": bool(power_diff[1] < power_diff[0] and power_diff[1]/abs(powers[-1]) <= .05),
        "continuum_budget_converges": all(abs(b["final_budget_defect"]) < abs(a["final_budget_defect"])
            for a, b in zip(temporal, temporal[1:])) and temporal[-1]["final_budget_over_actual_change"] <= .05,
        "weak_BE_work": all(r["max_weak_work_defect"] <= 1e-12 for r in experiments),
        "contact_topology": all(len(r["geometry"]["final"]["crossings"]) == 2 for r in experiments)}
    result = {"status": "complete", "qualification_status": "passed" if all(gates.values()) else "failed",
        "gates": gates, "mode": mode_record, "temporal": temporal, "amplitude_checks": amplitude_checks,
        "observed_modal_orders": orders.tolist(), "integrated_D_successive_differences": power_diff.tolist(),
        "amplitude_BE_error_sequence": amp_errors, "provenance": provenance,
        "runtime_s": time.perf_counter()-started, "scope": "isolated CH only, not full nonlinear CHNS qualification"}
    write_json(output/"status.json", result)
    print(json.dumps({"stage": "low_mode", "gates": gates, "orders": orders.tolist()}), flush=True)
    return result


def planner_stage():
    from ..ch_timestep_planner import ReducedCHReference, make_plan, verify_full_sparse_schedule
    op, _, operator_record, reference_record, basis, _ = qualified_reference()
    low = read_json(RESULTS/"low_mode/status.json")
    if (low.get("qualification_status") != "passed" or
            low["mode"]["operator_fingerprint"] != op.fingerprint or
            low["mode"]["Krylov_sha256"] != reference_record["Krylov_sha256"]):
        raise ValueError("STOP: matching low-mode prerequisite has not passed")
    output = reserve("stiff_bump/planner")
    started = time.perf_counter()
    V = basis["V"]
    chemical = op.mass_solve(op.H @ V)
    reduced = ReducedCHReference(basis["projected_L"], V.T @ (op.H @ V),
        op.mobility*chemical.T @ (op.K @ chemical),
        float(basis["beta"])*np.eye(V.shape[1])[:, 0], V.T @ (op.M0 @ V))
    fingerprints = {"operator": op.fingerprint,
        "equilibrium": operator_record["provenance"]["equilibrium_fingerprint"],
        "perturbation": operator_record["provenance"]["perturbation_coefficient_sha256"],
        "Krylov": reference_record["Krylov_sha256"]}
    plan = make_plan(reduced, 1e-4, reference_record["rho_relevant"], fingerprints)
    # The reduced candidate is always preserved, including all rejected levels.
    write_json(output/"reduced_candidate.json", plan)
    if plan["qualification_status"] != "linear_reduced_pass":
        result = {"status": "complete", "qualification_status": "failed",
            "reason": "STOP: no linear reduced candidate passed", "plan": plan,
            "runtime_s": time.perf_counter()-started}
        write_json(output/"status.json", result)
        return result
    seconds_per_step = operator_record["CH_dominated_justification"]["cost_seconds_per_step_conservative"]
    plan["Krylov_dimension"] = reference_record["dimension"]
    plan["cost_estimate"] = {"seconds_per_step": seconds_per_step, "level0_steps": plan["total_steps"],
        "level0_hours": plan["total_steps"]*seconds_per_step/3600,
        "level0_half_quarter_hours": 7*plan["total_steps"]*seconds_per_step/3600,
        "level0_limit_hours": 3., "series_limit_hours": 10.,
        "source": "conservative maximum of measured STEP3A.2 seconds per accepted interval"}
    cost_ok = plan["cost_estimate"]["level0_hours"] <= 3 and plan["cost_estimate"]["level0_half_quarter_hours"] <= 10
    plan["cost_gate"] = cost_ok
    sparse_ok = False
    if cost_ok:
        op.resolvent_backend = "petsc_mumps"
        def monitor(row):
            with (output/"full_sparse_progress.jsonl").open("a") as stream:
                stream.write(json.dumps(row)+"\n")
            print("Full sparse BE plan verification "+json.dumps(row), flush=True)
        verification = verify_full_sparse_schedule(op, reduced, V, plan["blocks"], monitor)
        write_json(output/"full_sparse_verification.json", verification)
        sparse_ok = (verification["endpoint_relative_L2_error"] <= 1e-3 and
            verification["max_checkpoint_relative_L2_error"] <= 1e-3 and
            verification["max_quadratic_energy_relative_error"] <= 1e-3 and
            verification["integrated_D2_relative_error"] <= 1e-2 and verification["max_mass_domain"] <= 1e-11)
        plan["full_sparse_verification"] = {k: v for k, v in verification.items() if k not in ("history", "blocks")}
        plan["full_sparse_verification_required"] = False
    else:
        plan["full_sparse_verification"] = {"status": "not_run", "reason": "predeclared production cost guard"}
    plan["nonlinear_execution_authorized"] = bool(cost_ok and sparse_ok)
    plan["qualification_status"] = "passed" if plan["nonlinear_execution_authorized"] else "failed"
    plan["plan_sha256"] = json_hash({k: v for k, v in plan.items() if k != "plan_sha256"})
    write_json(output/"timestep_plan.json", plan)
    result = {"status": "complete", "qualification_status": plan["qualification_status"],
        "reason": "All plan prerequisites passed" if plan["nonlinear_execution_authorized"] else
            ("STOP: predeclared production cost guard exceeded" if not cost_ok else "STOP: full sparse schedule failed"),
        "gates": {"reduced_BE": True, "cost": cost_ok, "full_sparse_BE": sparse_ok},
        "timestep_plan_sha256": hashlib.sha256((output/"timestep_plan.json").read_bytes()).hexdigest(),
        "semantic_plan_sha256": plan["plan_sha256"], "input_fingerprints": fingerprints,
        "dt0": plan["dt0"], "total_steps": plan["total_steps"], "cost_estimate": plan["cost_estimate"],
        "predicted_errors": plan["predicted_errors"], "runtime_s": time.perf_counter()-started}
    write_json(output/"status.json", result)
    print(json.dumps(result), flush=True)
    return result
