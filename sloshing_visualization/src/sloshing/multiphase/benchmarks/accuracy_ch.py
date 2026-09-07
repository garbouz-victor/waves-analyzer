"""STEP3A6: finite independent one-step accuracy experiments, no pilot/series."""
import copy
import gc
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile
import time
import xml.etree.ElementTree as ET
import numpy as np
from ..linearized_ch import array_hash, json_hash
from ..nonlinear_accuracy import (NonlinearControls, actual_controls, apply_controls,
    finite_arrays, synchronized_solution, check_candidate_material, publish_candidate,
    write_snapshot, read_snapshot)
from ..phase_rate import VERSION, ABSOLUTE_VERSION, PhaseMetric, IsolatedCHPhaseRate
from ..performance import Performance, petsc_events
from . import spectral_ch as historical
from .phase_rate_ch import HASHES, TINY_DT, MODERATE_DT, INITIAL_D, NewtonTrace, assemble

ROOT = Path("validation_results/step3a6")
P0 = Path("validation_results/step3a5")
CONTROLS = {"P0": NonlinearControls(), "P1": NonlinearControls(1e-12, 1e-12),
            "P2": NonlinearControls(3e-13, 3e-13)}
POLICY = {"moderate_dt": MODERATE_DT, "tiny_dt": TINY_DT, "delta": 1e-3,
    "controls": {k: v.record() for k, v in CONTROLS.items()}, "hashes": HASHES,
    "relative_comparison": 1e-10, "absolute_energy_mass_comparison": 1e-12,
    "newton_estimate_fraction": .25, "mass_domain": 1e-10, "certified_cells": 8.,
    "energy_growth": 1e-9, "weak_work": 1e-12, "linear_endpoint": 1e-10,
    "linear_increment": 1e-6, "linear_mass": 1e-11, "nonlinear_increment": 1e-3,
    "contrast_ratio": 10., "wall_budget_s": 1800.,
    "moderate_runs": ["P1_absolute", "P1_rate_semidiscrete", "P2_absolute", "P2_rate_semidiscrete", "P2_rate_zero"],
    "target_runs": ["rate_semidiscrete", "rate_zero", "absolute_control"],
    "pilot_authorized": False, "full_series_authorized": False, "full_CHNS_transfer": "not_qualified"}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def freeze():
    historical_policy = read(P0/"policy/policy.json")["policy"]
    if historical_policy["required_hashes"] != HASHES:
        raise ValueError("Frozen input policies disagree")
    checks = {"moderate_phi_mu_D_relative": POLICY["relative_comparison"],
        "moderate_energy_mass_absolute": POLICY["absolute_energy_mass_comparison"],
        "nonlinear_increment_relative": POLICY["nonlinear_increment"]}
    if any(historical_policy[k] != v for k, v in checks.items()):
        raise ValueError("Existing scientific thresholds changed")
    folder = ROOT/"policy"; folder.mkdir(parents=True, exist_ok=True)
    record = {"policy": POLICY, "policy_sha256": json_hash(POLICY),
        "design_sha256": hashlib.sha256(Path("STEP3A6_DESIGN.md").read_bytes()).hexdigest()}
    target = folder/"policy.json"
    if target.exists() and read(target) != record:
        raise ValueError("Frozen STEP3A6 policy/design changed")
    if not target.exists():
        write(target, record)
    return record


def budget_check(elapsed=0.):
    used = sum(read(p)["wall_s"] for p in ROOT.rglob("timing.json") if "tests" not in p.parts)
    if used+elapsed >= POLICY["wall_budget_s"]:
        raise RuntimeError("STOP: fixed 1800-second scientific budget")
    return used


def reserve(name):
    freeze(); budget_check()
    output = ROOT/name
    if ROOT.resolve() not in output.resolve().parents:
        raise ValueError("STEP3A6 output only")
    output.mkdir(parents=True, exist_ok=False)
    from ..provenance import source_hash
    from petsc4py import PETSc
    import dolfinx, scipy, ufl
    source = Path(__file__).resolve().parents[1]
    with tarfile.open(output/"source.tar.gz", "x:gz") as archive:
        for path in sorted(source.rglob("*.py")):
            archive.add(path, arcname=str(path.relative_to(source)))
    header = Path("/usr/local/petsc/linux-gnu-real64-32/include/dmumps_c.h").read_text()
    provenance = {"implementation": {"git_commit": os.environ.get("STEP3_GIT_COMMIT", "unavailable"),
        "source_sha256": source_hash(), "source_archive_sha256": hashlib.sha256((output/"source.tar.gz").read_bytes()).hexdigest(),
        "docker_image": os.environ.get("STEP3_DOCKER_IMAGE_DIGEST", "unavailable"),
        "rate_version": VERSION, "absolute_version": ABSOLUTE_VERSION},
        "historical_input_provenance": read(historical.RESULTS/"linear_operator/status.json")["provenance"],
        "policy_sha256": json_hash(POLICY), "required_hashes": HASHES,
        "stack": {"dolfinx": dolfinx.__version__, "PETSc": PETSc.Sys.getVersion(),
            "MUMPS": re.search(r'#define MUMPS_VERSION\s+"([^"]+)"', header).group(1),
            "scipy": scipy.__version__, "ufl": ufl.__version__, "scalar_dtype": str(np.dtype(PETSc.ScalarType)),
            "precision_bits": np.finfo(PETSc.ScalarType).bits}}
    write(output/"provenance.json", provenance)
    return output, provenance


def physical_checks(old, new, work, arrays):
    return {"finite": bool(finite_arrays(*arrays.values(), list(old.values()), list(new.values()),
                                        [v for v in work.values() if isinstance(v, (float, int))])),
        "mass": abs(new["mass_error_relative_to_domain"]) <= POLICY["mass_domain"],
        "resolution": new["cells_across_transition_certified_min"] >= POLICY["certified_cells"],
        "topology": new["wall_crossings"] == 2,
        "energy": new["E_total"]-old["E_total"] <= POLICY["energy_growth"],
        "BE_work": abs(work["weak_work_defect"]) <= POLICY["weak_work"],
        "decomposition": abs(work["decomposition_roundoff"]) <= POLICY["weak_work"]}


def comparison_gates(errors, absolute_physical, rate_physical):
    cross = {k: bool(np.isfinite(v) and v <= (POLICY["absolute_energy_mass_comparison"]
        if k in ("energy_absolute", "mass_absolute") else POLICY["relative_comparison"])) for k, v in errors.items()}
    return {"absolute_physical_checks": absolute_physical, "rate_physical_checks": rate_physical,
        "cross_formulation_checks": cross,
        "passed": all(absolute_physical.values()) and all(rate_physical.values()) and all(cross.values())}


def may_continue_accuracy(result):
    """P1 cross FAIL is not a physical FAIL; a failed primary never authorizes tiny."""
    return result["nonlinear_solver_status"] == "converged" and result["physical_checks"] == "passed"


def sample(observer, diagnostics):
    from ..diagnostics import interface_resolution
    from ..ch_validation_diagnostics import WallTrace
    row = diagnostics.measure(); resolution = interface_resolution(observer)
    row.update(cells_across_transition_certified_min=resolution.get("cells_across_transition_certified_min", 0.),
        wall_crossings=len(WallTrace(observer).crossings()))
    return row, resolution


def residual_audit(solver, arrays, metric, space, *, include_vectors=False):
    """Private archive reconstruction; never changes accepted state or runs Newton."""
    import ufl
    from dolfinx import fem
    from ..equilibrium import _measures, chemical_stationary_form
    pmap = np.asarray(space.sub(0).collapse()[1]); mmap = np.asarray(space.sub(1).collapse()[1])
    old, physical = fem.Function(space), fem.Function(space)
    old.x.array[pmap] = arrays["phi_old"]; old.x.array[mmap] = arrays["mu_old"]
    physical.x.array[pmap] = arrays["phi_new"]; physical.x.array[mmap] = arrays["mu_new"]
    old.x.scatter_forward(); physical.x.scatter_forward()
    phi, mu = ufl.split(physical); phi0, _ = ufl.split(old)
    test, chi = ufl.TestFunctions(space)
    dx, _ = _measures(solver); dt = float(arrays["dt"]); mob = solver.config.mobility
    vectors = {"literal_absolute_UFL": assemble(((phi-phi0)/dt*test+mob*ufl.inner(ufl.grad(mu),ufl.grad(test)))*dx)[pmap],
        "coefficient_difference_sparse": metric.M @ ((arrays["phi_new"]-arrays["phi_old"])/dt)+mob*(metric.K @ arrays["mu_new"]),
        "chemical_reconstructed": assemble(mu*chi*dx-chemical_stationary_form(solver,phi,0.,chi))[mmap]}
    if "phase_rate" in arrays:
        unknown = fem.Function(space)
        unknown.x.array[pmap] = arrays["phase_rate"]; unknown.x.array[mmap] = arrays["mu_new"]
        unknown.x.scatter_forward()
        r, chemical = ufl.split(unknown)
        old_scalar = fem.Function(space.sub(0).collapse()[0])
        old_scalar.x.array[:] = arrays["phi_old"]; old_scalar.x.scatter_forward()
        # Constant dt matches the actual solver expression, not a post-hoc rate.
        from petsc4py import PETSc
        dt_coefficient = fem.Constant(solver.mesh, PETSc.ScalarType(dt))
        vectors.update(retained_rate_UFL=assemble((r*test+mob*ufl.inner(ufl.grad(chemical),ufl.grad(test)))*dx)[pmap],
            retained_rate_sparse=metric.M @ arrays["phase_rate"]+mob*(metric.K @ arrays["mu_new"]),
            chemical_rate_expression=assemble(chemical*chi*dx-chemical_stationary_form(solver,
                old_scalar+dt_coefficient*r,0.,chi))[mmap])
    norms = {k: {"algebraic_norm": float(np.linalg.norm(v)), "Riesz_L2": metric.riesz_norm(v)} for k,v in vectors.items()}
    return (norms, vectors) if include_vectors else norms


def newton_estimate(s, engine, arrays, metric, perf):
    """ONE sparse J e=-R, never publish e; no rigorous-error-bound claim."""
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from scipy.sparse import csr_matrix
    from ..linearized_ch import PETScSparseFactor
    from ..equilibrium import _measures, chemical_stationary_form
    space = engine.state.function_space
    pmap = np.asarray(space.sub(0).collapse()[1]); mmap = np.asarray(space.sub(1).collapse()[1])
    with perf.measure("assembly"):
        residual = assemble(engine.F)
        matrix = fp.assemble_matrix(fem.form(ufl.derivative(engine.F,engine.state,ufl.TrialFunction(space))))
        matrix.assemble(); ptr,col,val = matrix.getValuesCSR()
        J = csr_matrix((val.copy(),col.copy(),ptr.copy()),shape=matrix.getSize()); matrix.destroy()
    with perf.measure("diagnostic_Newton_factor"):
        factor = PETScSparseFactor(J)
    with perf.measure("diagnostic_Newton_solve"):
        e = factor.solve(-residual); reason = factor.ksp.getConvergedReason()
    factor.close()
    ep = e[pmap]*float(arrays["dt"]) if "phase_rate" in arrays else e[pmap]
    emu = e[mmap]; mu = arrays["mu_new"]; mob = s.config.mobility
    D = mob*float(mu @ (metric.K @ mu))
    dD = 2*mob*float(mu @ (metric.K @ emu)); quadratic = mob*float(emu @ (metric.K @ emu))
    Q = space.sub(0).collapse()[0]; p, dp = fem.Function(Q), fem.Function(Q)
    p.x.array[:] = arrays["phi_new"]; dp.x.array[:] = ep
    p.x.scatter_forward(); dp.x.scatter_forward()
    energy = float(fem.assemble_scalar(fem.form(chemical_stationary_form(s,p,0.,dp))))
    errors = {"phi_relative": metric.norm(ep)/metric.norm(arrays["phi_new"]),
        "mu_relative": metric.norm(emu)/metric.norm(mu), "D_relative": (abs(dD)+abs(quadratic))/abs(D),
        "energy_absolute": abs(energy), "mass_absolute": abs(float(metric.mass @ ep))}
    gates = {k: bool(np.isfinite(v) and v <= .25*(POLICY["absolute_energy_mass_comparison"]
             if k.endswith("absolute") else POLICY["relative_comparison"])) for k,v in errors.items()}
    return {"estimated_errors": errors, "engineering_checks": gates, "D_linear_correction": dD,
        "D_quadratic_term": quadratic, "linear_solve_reason": reason,
        "linear_solve_relative_residual": float(np.linalg.norm(J @ e+residual)/max(np.linalg.norm(residual),1e-300)),
        "finite": bool(finite_arrays(e)), "published": False, "rigorous_bound": False}


def run_one(name, formulation, guess, level, dt, *, estimate=False):
    from dolfinx import fem
    from petsc4py import PETSc
    from ..isolated_ch import IsolatedCH
    from ..diagnostics import Diagnostics
    from ..be_work_diagnostics import BEWorkDiagnostics
    from ..equilibrium import equilibrium_config_fingerprint
    perf = Performance(); output, provenance = reserve(name); PETSc.Log.begin()
    engine, trace, evidence, accepted_before = None, None, {}, None
    result = {"execution_status": "failed", "nonlinear_solver_status": "not_run", "physical_checks": "not_run",
        "accepted": False, "formulation": formulation, "initial_guess": guess, "accuracy_level": level, "dt": dt}
    try:
        with perf.measure("initialization_JIT"):
            s, eq, psi, initial, input_provenance = historical.scientific_input()
            field = fem.Function(eq.function_space); field.x.array[:] = initial; field.x.scatter_forward()
            engine = IsolatedCH(s,field) if formulation == "absolute" else IsolatedCHPhaseRate(s,field,initial_guess=guess)
            engine.dt.value = dt; engine.performance = perf
            controls = CONTROLS[level]
            apply_controls(engine.problem.solver, controls, "step3a6_"+name.replace("/","_")+"_")
            trace = NewtonTrace(engine,output)
            before_controls = actual_controls(engine.problem.solver)
            if before_controls != controls.record():
                raise ValueError("Actual controls differ before solve")
            old_phi, old_mu = s.state.sub(2).collapse().x.array.copy(),s.state.sub(3).collapse().x.array.copy()
            metric = engine.metric if formulation == "rate" else PhaseMetric(s,eq.function_space)
            pmap = np.asarray(engine.state.function_space.sub(0).collapse()[1])
            mmap = np.asarray(engine.state.function_space.sub(1).collapse()[1])
            guess_phi = (old_phi+dt*engine.state.x.array[pmap] if formulation == "rate" else engine.state.x.array[pmap].copy())
            initial_checks = {"initial_phi": array_hash(old_phi)==HASHES["initial_phi_sha256"],
                "psi": array_hash(psi)==HASHES["perturbation_coefficient_sha256"],
                "equilibrium": input_provenance["equilibrium_fingerprint"]==HASHES["equilibrium_fingerprint"]}
            plan = read(Path("validation_results/step3a4/schedule_optimization/optimized_plan.json"))
            initial_checks["plan"] = (plan["plan_sha256"]==HASHES["optimized_plan_sha256"] and
                json_hash({k:v for k,v in plan.items() if k!="plan_sha256"})==plan["plan_sha256"])
            if not all(initial_checks.values()):
                raise ValueError("Initial physical/hash mismatch")
            provenance.update(physical_config=s.config.as_dict(), physical_config_fingerprint=equilibrium_config_fingerprint(s.config),
                mesh_fingerprint=input_provenance["mesh_fingerprint"], nonlinear_solver_controls=before_controls,
                actual_initial_physical_guess_sha256=array_hash(guess_phi), old_physical_phi_sha256=array_hash(old_phi),
                mu_old_sha256=array_hash(old_mu), dt=dt, initial_guess_policy=guess)
            initial_arrays=dict(phi_old=old_phi,mu_old=old_mu,phi_guess=guess_phi,mixed_guess=engine.state.x.array.copy())
            provenance["initial_array_sha256"]={k:array_hash(v) for k,v in initial_arrays.items()}
            write(output/"provenance.json",provenance)
            np.savez_compressed(output/"initial.npz",**initial_arrays)
            old_diag = Diagnostics(s)
        with perf.measure("diagnostics"):
            old, old_resolution = sample(s,old_diag)
        initial_checks.update(D0=abs(old["CH_dissipation"]/INITIAL_D-1)<=1e-10,
            resolution=old["cells_across_transition_certified_min"]>=8, topology=old["wall_crossings"]==2)
        write(output/"initial.json",{"state":old,"checks":initial_checks,"resolution":old_resolution,
            "guess_scale": engine.rate_statistics(dt) if formulation=="rate" else None})
        if not all(initial_checks.values()):
            raise RuntimeError("Initial scientific gate failed")
        accepted_before = {k:getattr(s,k).x.array.copy() for k in ("state","old","older")}
        accepted_attrs = {k:getattr(s,k) for k in ("time","step_number","current_dt","current_phase")}
        accepted_engine = {k:copy.deepcopy(getattr(engine,k)) for k in ("history","last_accepted_rate","last_snapshot") if hasattr(engine,k)}
        def validate(candidate, arrays, candidate_dt, end_time):
            observer = copy.copy(s); observer.state=candidate; observer.older=s.state
            observer.time=end_time; observer.step_number=s.step_number+1; observer.current_dt=candidate_dt
            with perf.measure("diagnostics"):
                diagnostic = Diagnostics(observer)
                diagnostic.initial=copy.deepcopy(old_diag.initial); diagnostic.previous=copy.deepcopy(old_diag.previous)
                diagnostic.cumulative=old_diag.cumulative.copy()
                new,resolution = sample(observer,diagnostic)
                work = BEWorkDiagnostics(observer).measure(candidate_dt,old,new)
                checks = physical_checks(old,new,work,arrays)
                checks["material"]=True  # same existing checker already passed on this private candidate
            evidence.update(old=old,new=new,resolution=resolution,BE_work=work,checks=checks)
            if not all(checks.values()):
                raise RuntimeError("Physical candidate gate failed; no publication")
        if formulation == "rate":
            log=engine.step(dt,candidate_validator=validate)
            arrays={k:v.copy() for k,v in engine.last_snapshot.items()}
        else:
            started=time.perf_counter()
            with perf.measure("nonlinear_SNES"):
                engine.problem.solve()
            synchronized_solution(engine.problem,engine.state)
            arrays={"phi_old":old_phi,"mu_old":old_mu,"phi_new":engine.state.x.array[pmap].copy(),
                "mu_new":engine.state.x.array[mmap].copy(),"dt":np.array(dt),
                "mixed_Function":engine.state.x.array.copy(),"PETSc_solution":engine.problem.x.array.copy(),
                "rate_map":pmap,"mu_map":mmap}
            arrays["increment_from_physical_states"]=arrays["phi_new"]-old_phi
            candidate=fem.Function(s.space); candidate.x.array[:]=s.state.x.array
            candidate.x.array[np.asarray(s.space.sub(2).collapse()[1])]=arrays["phi_new"]
            candidate.x.array[np.asarray(s.space.sub(3).collapse()[1])]=arrays["mu_new"]
            candidate.x.scatter_forward(); check_candidate_material(s,candidate)
            validate(candidate,arrays,dt,dt)
            publish_candidate(s,candidate,dt,dt,"isolated_absolute_accuracy_BE")
            log={"runtime_s":time.perf_counter()-started}
        arrays["mu_old"]=old_mu.copy()
        after_controls=actual_controls(engine.problem.solver)
        if after_controls != before_controls:
            raise RuntimeError("Controls changed during solve")
        with perf.measure("archival_IO"):
            # Archive is provisional until roundtrip and read-only checks finish.
            snapshot=write_snapshot(output,arrays,{"accepted":False,"execution_status":"pending_validation","dt":dt,
                "layout":"serial native P2/P2 mixed; scalar collapse maps stored", "formulation_version":VERSION if formulation=="rate" else ABSOLUTE_VERSION,
                "mesh_fingerprint":provenance["mesh_fingerprint"], "provenance":provenance})
            restored,metadata=read_snapshot(output)
        roundtrip=all(np.array_equal(restored[k],v) for k,v in arrays.items())
        if not roundtrip:
            raise RuntimeError("Retained archive roundtrip failed")
        with perf.measure("residual_reassembly"):
            audit,vectors=residual_audit(s,restored,metric,engine.state.function_space,include_vectors=True)
            raw_phase=assemble(engine.phase_form)[pmap]; raw_chemical=assemble(engine.chemical_form)[mmap]
            vectors.update(native_phase=raw_phase,native_chemical=raw_chemical)
            np.savez_compressed(output/"residual_vectors.npz",**vectors)
        raw={"phase_algebraic":float(np.linalg.norm(raw_phase)),"chemical_algebraic":float(np.linalg.norm(raw_chemical)),
            "phase_Riesz_L2":metric.riesz_norm(raw_phase),"chemical_Riesz_L2":metric.riesz_norm(raw_chemical)}
        error_estimate=newton_estimate(s,engine,arrays,metric,perf) if estimate else None
        snapshot.update(accepted=True,execution_status="complete")
        snapshot["snapshot_sha256"]=json_hash({k:v for k,v in snapshot.items() if k!="snapshot_sha256"})
        write(output/"snapshot.json",snapshot)
        result.update(execution_status="complete",nonlinear_solver_status="converged",physical_checks="passed",accepted=True,
            physical_check_details=evidence["checks"],evidence=evidence,SNES=trace.summary(), controls_before=before_controls,
            controls_after=after_controls, residuals=audit, actual_solution_residual=raw,
            newton_error_estimate=error_estimate, snapshot_sha256=snapshot["snapshot_sha256"],
            archive_roundtrip=roundtrip, retained_rate_sha256=snapshot["array_sha256"].get("phase_rate"),
            residual_vector_sha256={k:array_hash(v) for k,v in vectors.items()},
            reconstructed_phi_sha256=snapshot["array_sha256"]["phi_new"], log=log)
        write(output/"status.json",result)
        print(name+" "+json.dumps({"SNES":result["SNES"],"physical":result["physical_checks"],"estimate":error_estimate}),flush=True)
        return result
    except Exception as error:
        result.update(execution_status="failed",accepted=False,error=str(error), numerical_accuracy_floor_detected=None,
            accuracy_floor_note="A diverged solve alone does not prove a floating-point floor")
        if accepted_before is not None:
            for k,v in accepted_before.items():
                getattr(s,k).x.array[:]=v; getattr(s,k).x.scatter_forward()
            for k,v in accepted_attrs.items(): setattr(s,k,v)
            for k in ("history","last_accepted_rate","last_snapshot"):
                if k in accepted_engine: setattr(engine,k,accepted_engine[k])
                elif hasattr(engine,k): delattr(engine,k)
        if (output/"snapshot.json").exists():
            rejected_snapshot=read(output/"snapshot.json")
            rejected_snapshot.update(accepted=False,execution_status="failed")
            rejected_snapshot["snapshot_sha256"]=json_hash({k:v for k,v in rejected_snapshot.items() if k!="snapshot_sha256"})
            write(output/"snapshot.json",rejected_snapshot)
        if engine is not None:
            result["SNES"]=trace.summary() if trace else None
            reason=engine.problem.solver.getConvergedReason()
            result["nonlinear_solver_status"]="converged" if reason>0 else "diverged"
            result["physical_checks"]="failed" if reason>0 else "not_run"
            rejected={"last_evaluated_Function":engine.state.x.array.copy(),
                "internal_SNES_iterate":engine.problem.x.array.copy(),
                "accepted_physical_state":s.state.x.array.copy(),
                "phi_old":s.state.sub(2).collapse().x.array.copy(),
                "mu_old":s.state.sub(3).collapse().x.array.copy(),
                "phase_map":np.asarray(engine.state.function_space.sub(0).collapse()[1]),
                "mu_map":np.asarray(engine.state.function_space.sub(1).collapse()[1]),"dt":np.array(dt)}
            np.savez_compressed(output/"rejected_attempt.npz",**rejected)
            write(output/"rejected_attempt.json",{"accepted":False,"execution_status":"failed",
                "array_sha256":{k:array_hash(v) for k,v in rejected.items()},
                "note":"Function last residual trial, internal SNES iterate, and physical accepted state are DISTINCT; no rejected field published"})
        if evidence:
            result["physical_check_details"]=evidence.get("checks",{})
            result["evidence"]=evidence
        write(output/"status.json",result); print(name+" STOP "+str(error),flush=True)
        return result
    finally:
        write(output/"timing.json",{**perf.snapshot(),"PETSc":petsc_events()})
        gc.collect()


def compare_runs(absolute_path, rate_path, metric):
    a,b=read(absolute_path/"status.json"),read(rate_path/"status.json")
    aa,_=read_snapshot(absolute_path); bb,_=read_snapshot(rate_path)
    errors={k+"_relative":metric.norm(bb[k+"_new"]-aa[k+"_new"])/metric.norm(aa[k+"_new"]) for k in ("phi","mu")}
    x,y=a["evidence"]["new"],b["evidence"]["new"]
    errors.update(D_relative=abs(y["CH_dissipation"]/x["CH_dissipation"]-1),
        energy_absolute=abs(y["E_total"]-x["E_total"]),mass_absolute=abs(y["phase_mass"]-x["phase_mass"]))
    return {"errors":errors,**comparison_gates(errors,a["physical_check_details"],b["physical_check_details"])}


def moderate_accuracy():
    freeze()
    test_path=ROOT/"tests/cheap_docker_final.xml"
    test_root=ET.parse(test_path).getroot()
    required={"test_tiny_FE_residual_and_jacobian_equivalence_and_FD",
        "test_strict_tiny_FEM_moderate_comparison_and_matched_guess"}
    successful={case.attrib["name"] for case in test_root.iter("testcase") if len(case)==0}
    if not required <= successful or list(test_root.iter("failure")) or list(test_root.iter("error")):
        raise ValueError("STOP: required cheap pinned algebraic/Jacobian tests not passed")
    test_evidence={"path":str(test_path),"sha256":hashlib.sha256(test_path.read_bytes()).hexdigest(),"passed":True}
    folder=ROOT/"moderate_accuracy"; folder.mkdir(exist_ok=True)
    results={}
    for level in ("P1","P2"):
        for formulation,guess in (("absolute","phi_old"),("rate","semidiscrete")):
            name=level+"_"+("absolute" if formulation=="absolute" else "rate_semidiscrete")
            results[name]=run_one("moderate_accuracy/"+name,formulation,guess,level,MODERATE_DT,estimate=level=="P2")
            if not may_continue_accuracy(results[name]):
                write(folder/"comparison.json",{"moderate_accuracy":"inconclusive","stop_run":name,"runs":results})
                return
    results["P2_rate_zero"]=run_one("moderate_accuracy/P2_rate_zero","rate","zero","P2",MODERATE_DT,estimate=True)
    if not may_continue_accuracy(results["P2_rate_zero"]):
        write(folder/"comparison.json",{"moderate_accuracy":"inconclusive","stop_run":"P2_rate_zero","runs":results})
        return
    op,_,_=historical.load_operator()
    comparisons={level:compare_runs(folder/(level+"_absolute"),folder/(level+"_rate_semidiscrete"),op) for level in ("P1","P2")}
    matched=compare_runs(folder/"P2_absolute",folder/"P2_rate_zero",op)
    pa,pz=(read(folder/name/"provenance.json") for name in ("P2_absolute","P2_rate_zero"))
    matched_inputs=all(pa[key]==pz[key] for key in ("actual_initial_physical_guess_sha256",
        "old_physical_phi_sha256","mu_old_sha256","nonlinear_solver_controls","physical_config_fingerprint","dt"))
    matched["same_physical_initial_guess_and_inputs"]=matched_inputs
    matched["initial_residual_absolute_difference"]=abs(results["P2_absolute"]["SNES"]["initial_residual"]-
        results["P2_rate_zero"]["SNES"]["initial_residual"])
    accuracy={k:v["newton_error_estimate"] for k,v in results.items() if v["newton_error_estimate"] is not None}
    accuracy_ok=all(all(v["engineering_checks"].values()) and v["finite"] and v["linear_solve_reason"]>0 for v in accuracy.values())
    qualified=comparisons["P2"]["passed"] and matched["passed"] and matched_inputs and accuracy_ok and all(v["archive_roundtrip"] for v in results.values())
    changes={}
    for formulation,old_name in (("absolute","absolute"),("rate_semidiscrete","phase_rate")):
        with np.load(P0/"moderate_dt"/old_name/"physical.npz") as z:
            previous={k+"_new":z[k].copy() for k in ("phi","mu")}
        previous_D=read(P0/"moderate_dt/status.json")["runs"][old_name]["new"]["CH_dissipation"]
        for level in ("P1","P2"):
            current,_=read_snapshot(folder/(level+"_"+formulation))
            D=results[level+"_"+formulation]["evidence"]["new"]["CH_dissipation"]
            changes[level+"_"+formulation]={k+"_relative":op.norm(current[k+"_new"]-previous[k+"_new"])/op.norm(current[k+"_new"]) for k in ("phi","mu")}
            changes[level+"_"+formulation]["D_relative"]=abs(previous_D/D-1)
            previous,previous_D=current,D
    write(folder/"newton_error_estimates.json",accuracy)
    record={"moderate_accuracy":"qualified" if qualified else "failed", "comparisons":comparisons,
        "matched_guess_control":matched,"nonlinear_accuracy_checks":{"engineering_estimates":accuracy_ok},
        "within_formulation_tighter_changes":changes,"algebraic_equivalence_test_evidence":test_evidence,
        "P0_historical":read(P0/"moderate_dt/status.json"),"runs":results}
    write(folder/"comparison.json",record); print("MODERATE "+record["moderate_accuracy"],flush=True)


def require_moderate():
    if read(ROOT/"moderate_accuracy/comparison.json")["moderate_accuracy"]!="qualified":
        raise ValueError("STOP: moderate accuracy not qualified")


def linear_tiny():
    require_moderate()
    from scipy.sparse import bmat
    from ..linearized_ch import PETScSparseFactor
    from petsc4py import PETSc
    perf=Performance(); output,provenance=reserve("linear_tiny"); PETSc.Log.begin()
    try:
        with perf.measure("initialization"):
            op,q0,record=historical.load_operator()
            op.resolvent_backend,op.performance="petsc_mumps",perf
            with np.load(historical.RESULTS/"linear_operator/inputs.npz") as z:
                psi=z["psi"].copy()
            if array_hash(psi)!=HASHES["perturbation_coefficient_sha256"] or not np.array_equal(q0,op.project(psi)):
                raise ValueError("Historical q0 convention/hash mismatch")
        qnew=op.be_step(q0,TINY_DT)
        with perf.measure("assembly"):
            matrix=bmat([[op.M0,op.mobility*op.K],[-TINY_DT*op.H,op.M0]],format="csr")
            rhs=np.concatenate((np.zeros(op.size),op.H @ q0))
        with perf.measure("factorization_setup"):
            factor=PETScSparseFactor(matrix)
        with perf.measure("linear_solve"):
            solution=factor.solve(rhs); reason=factor.ksp.getConvergedReason()
        factor.close(); op.clear_be_cache()
        rate,mu=solution[:op.size],solution[op.size:]
        qrate=q0+TINY_DT*rate
        errors={"endpoint_relative":op.norm(qrate-qnew)/op.norm(qnew),
            "increment_relative":op.norm(TINY_DT*rate-(qnew-q0))/op.norm(qnew-q0),
            "mass_domain":max(abs(float(op.mass @ (x-q0)))/op.area for x in (qnew,qrate))}
        checks={k:bool(np.isfinite(v) and v<=POLICY["linear_"+k.split("_")[0]]) for k,v in errors.items()}
        checks["finite"]=finite_arrays(solution,qnew); checks["KSP_converged"]=reason>0
        arrays=dict(q0=q0,q_new_old_linear=qnew,q_new_rate_linear=qrate,rate_linear=rate,mu_linear=mu,
            increment_from_rate=TINY_DT*rate,increment_old_physical=qnew-q0,
            increment_rate_physical=qrate-q0,dt=np.array(TINY_DT))
        np.savez_compressed(output/"linear_step.npz",**arrays)
        provenance.update(operator_fingerprint=op.fingerprint,dt=TINY_DT,
            q0_convention="mass-projected historical psi WITHOUT delta",nonlinear_delta=POLICY["delta"])
        write(output/"provenance.json",provenance)
        status={"execution_status":"complete","linear_comparison":"passed" if all(checks.values()) else "failed",
            "benchmark":"full sparse LINEAR one-step only", "errors":errors,"checks":checks,
            "array_sha256":{k:array_hash(v) for k,v in arrays.items()},"KSP_reason":reason,
            "sparse_system_relative_residual":float(np.linalg.norm(matrix @ solution-rhs)/np.linalg.norm(rhs))}
        write(output/"status.json",status); print("LINEAR TINY "+json.dumps(status),flush=True)
        return status
    finally:
        write(output/"timing.json",{**perf.snapshot(),"PETSc":petsc_events()})


def target_diagnostics(name):
    """Archive-only arithmetic; all FEM residual vectors were saved by run_one."""
    from decimal import Decimal, localcontext
    budget_check(); perf=Performance()
    output=ROOT/"target_tiny"/name
    arrays,metadata=read_snapshot(output); status=read(output/"status.json")
    if not metadata["accepted"] or not status["accepted"]:
        raise ValueError("No proof from a rejected state")
    op,_,_=historical.load_operator()
    linear_status=read(ROOT/"linear_tiny/status.json")
    with np.load(ROOT/"linear_tiny/linear_step.npz") as z:
        linear={k:z[k].copy() for k in z.files}
    if {k:array_hash(v) for k,v in linear.items()}!=linear_status["array_sha256"]:
        raise ValueError("Linear archive hash mismatch")
    retained=arrays["increment_from_rate"]; stored=arrays["increment_from_physical_states"]
    discrepancy=stored-retained
    rounding={"increment_rate_L2":op.norm(retained),"increment_stored_L2":op.norm(stored),
        "discrepancy_L2":op.norm(discrepancy),"discrepancy_max":float(np.max(abs(discrepancy))),
        "discrepancy_relative":op.norm(discrepancy)/op.norm(retained),
        "max_coefficient_ULP_over_dt":float(np.max(abs(np.spacing(arrays["phi_new"]))))/TINY_DT,
        "bitwise_reconstruction":bool(np.array_equal(arrays["phi_old"]+TINY_DT*arrays["phase_rate"],arrays["phi_new"]))}
    selected=np.argsort(-abs(discrepancy),kind="stable")[:5]
    decimal=[]
    with localcontext() as context:
        context.prec=80
        for i in selected:
            exact=Decimal.from_float(float(arrays["phi_old"][i]))+Decimal.from_float(TINY_DT)*Decimal.from_float(float(arrays["phase_rate"][i]))
            error=Decimal.from_float(float(arrays["phi_new"][i]))-exact
            decimal.append({"dof":int(i),"exact_binary_input_sum_decimal":str(exact),
                "rounded_minus_exact":str(error),"ULP":float(abs(np.spacing(arrays["phi_new"][i])))})
    rounding["selected_Decimal_audit"]=decimal
    delta=POLICY["delta"]
    linear_inc=delta*linear["increment_old_physical"]
    linear_retained=delta*linear["increment_from_rate"]
    mu_reference=arrays["mu_old"]+delta*(linear["mu_linear"]-op.chemical(linear["q0"]))
    comparison={"increment_absolute_L2":op.norm(retained-linear_inc),
        "increment_relative_L2":op.norm(retained-linear_inc)/op.norm(linear_inc),
        "retained_linear_increment_relative_L2":op.norm(retained-linear_retained)/op.norm(linear_retained),
        "endpoint_relative_L2":op.norm(arrays["phi_new"]-(arrays["phi_old"]+linear_inc))/op.norm(arrays["phi_new"]),
        "chemical_update_relative_L2_diagnostic":op.norm(arrays["mu_new"]-mu_reference)/max(op.norm(mu_reference-arrays["mu_old"]),1e-300),
        "q0_convention":"psi without delta; delta multiplied exactly once",
        "threshold":POLICY["nonlinear_increment"]}
    comparison["passed"]=comparison["increment_relative_L2"]<=POLICY["nonlinear_increment"]
    old,new,work=(status["evidence"][k] for k in ("old","new","BE_work"))
    physics={"mass_domain":new["mass_error_relative_to_domain"],
        "certified_cells":new["cells_across_transition_certified_min"],
        "phi_min":new["phi_min_dof"],"phi_max":new["phi_max_dof"],"mu_mean":new["mu_mean_Pa"],"mu_std":new["mu_std_Pa"],
        "F_old":old["E_total"],"F_new":new["E_total"],
        "DeltaF_total_subtraction":new["E_total"]-old["E_total"],
        "DeltaF_existing_integrand_difference":work["Delta_E"],
        "DeltaF_existing_work_minus_remainder":work["energy_work"]-work["B_BE"],
        "weak_work_relative_to_interval_activity_diagnostic":abs(work["weak_work_defect"])/max(abs(work["energy_work"]-work["B_BE"]),work["trapezoid_integral"],1e-30),
        "BE_work":work,"rate_L2":op.norm(arrays["phase_rate"]),
        "rate_max":float(np.max(abs(arrays["phase_rate"]))),"dt_rate_max":float(np.max(abs(retained)))}
    rejected_path=Path("validation_results/step3a4/performance/isolated_pilot/failed_iterate.npz")
    historic_status=read(Path("validation_results/step3a4/performance/failed_interval_audit/status.json"))
    if historic_status["dt"]!=TINY_DT:
        raise ValueError("Historical rejected dt mismatch")
    with np.load(rejected_path) as z:
        # This archive is a Function residual-trial, NOT an accepted solution.
        if array_hash(z["initial_phi"])!=HASHES["initial_phi_sha256"]:
            raise ValueError("Historical rejected input mismatch")
        rejected={"label":"STEP3A4 REJECTED last Function, not a solution", "dt":TINY_DT,
            "phi_difference_L2":op.norm(arrays["phi_new"]-z["failed_phi"]),
            "mu_difference_L2":op.norm(arrays["mu_new"]-z["failed_mu"]),
            "phi_max_coefficient_difference":float(np.max(abs(arrays["phi_new"]-z["failed_phi"]))),
            "D_new_difference":new["CH_dissipation"]-historic_status["rejected_state"]["CH_dissipation"],
            "DeltaF_total_difference":physics["DeltaF_total_subtraction"]-
                (historic_status["rejected_state"]["E_total"]-historic_status["initial"]["E_total"]),
            "DeltaF_integrand_difference":work["Delta_E"]-historic_status["rejected_BE_work_diagnostic"]["Delta_E"],
            "historical_audit":historic_status}
    result={"residuals":status["residuals"],"reconstruction":rounding,"linear_comparison":comparison,
        "physics":physics,"historical_rejected_comparison":rejected}
    write(output/"proof_diagnostics.json",result)
    (output/"read_only").mkdir()
    write(output/"read_only/timing.json",perf.snapshot())
    return result


def target_tiny():
    require_moderate()
    if read(ROOT/"linear_tiny/status.json")["linear_comparison"]!="passed":
        raise ValueError("STOP: full-sparse tiny comparison not passed")
    folder=ROOT/"target_tiny"; folder.mkdir(exist_ok=True)
    primary=run_one("target_tiny/rate_semidiscrete","rate","semidiscrete","P0",TINY_DT)
    result={"execution_status":primary["execution_status"],"target_tiny_step":"rejected",
        "cancellation_hypothesis":"inconclusive","runs":{"rate_semidiscrete":primary},
        "pilot_authorized":False,"full_series_authorized":False,"full_CHNS_transfer":"not_qualified"}
    if not may_continue_accuracy(primary):
        write(folder/"comparison.json",result); return result
    result["target_tiny_step"]="accepted"
    proof=target_diagnostics("rate_semidiscrete"); result["primary_proof"]=proof
    if not proof["linear_comparison"]["passed"]:
        result["stop_reason"]="Nonlinear/sparse linear increment comparison failed"
        write(folder/"comparison.json",result); return result
    zero=run_one("target_tiny/rate_zero","rate","zero","P0",TINY_DT)
    result["runs"]["rate_zero"]=zero
    if zero["nonlinear_solver_status"]=="converged" and zero["physical_checks"]!="passed":
        result["stop_reason"]="Matched-guess physical failure blocks further PDE work"
        write(folder/"comparison.json",result); return result
    if zero["accepted"]:
        result["zero_proof"]=target_diagnostics("rate_zero")
        op,_,_=historical.load_operator()
        result["rate_guess_comparison"]=compare_runs(folder/"rate_semidiscrete",folder/"rate_zero",op)
    absolute=run_one("target_tiny/absolute_control","absolute","phi_old","P0",TINY_DT)
    result["runs"]["absolute_control"]=absolute
    pp,pz,pa=(read(folder/k/"provenance.json") for k in ("rate_semidiscrete","rate_zero","absolute_control"))
    inputs=all(pp[key]==pz[key]==pa[key] for key in ("old_physical_phi_sha256","mu_old_sha256","nonlinear_solver_controls","physical_config_fingerprint","dt"))
    matched=(pz["actual_initial_physical_guess_sha256"]==pa["actual_initial_physical_guess_sha256"])
    rate_norm=primary["residuals"]["retained_rate_UFL"]["algebraic_norm"]
    absolute_norm=primary["residuals"]["literal_absolute_UFL"]["algebraic_norm"]
    contrast=absolute_norm>POLICY["contrast_ratio"]*rate_norm and absolute_norm>primary["SNES"]["effective_target"]
    supported=(inputs and matched and zero["accepted"] and result["zero_proof"]["linear_comparison"]["passed"] and
        result["rate_guess_comparison"]["passed"] and
        absolute["nonlinear_solver_status"]=="diverged" and contrast and proof["reconstruction"]["bitwise_reconstruction"])
    hypothesis="supported" if supported else ("not_reproduced" if absolute["accepted"] else "inconclusive")
    result.update(cancellation_hypothesis=hypothesis,same_inputs=inputs,matched_zero_absolute_initial_guess=matched,
        residual_contrast=contrast,literal_to_retained_ratio=absolute_norm/max(rate_norm,1e-300))
    write(folder/"comparison.json",result); print("TARGET TINY "+hypothesis,flush=True)
    return result


def analyze():
    """Read-only result aggregation; no solver or downstream authorization."""
    import csv
    moderate_path=ROOT/"moderate_accuracy/comparison.json"
    target_path=ROOT/"target_tiny/comparison.json"
    moderate=read(moderate_path) if moderate_path.exists() else {"moderate_accuracy":"not_run"}
    target=read(target_path) if target_path.exists() else {"target_tiny_step":"not_run","cancellation_hypothesis":"inconclusive"}
    if moderate["moderate_accuracy"]!="qualified":
        verdict="MODERATE NONLINEAR ACCURACY NOT QUALIFIED; TARGET TINY STEP NOT RUN"
    elif target["target_tiny_step"]!="accepted":
        verdict="MODERATE FORMULATION EQUIVALENCE QUALIFIED; TINY-dt BLOCKER REMAINS"
    elif target["cancellation_hypothesis"]=="supported":
        verdict="TARGET TINY-dt ISOLATED CH BE STEP SOLVED; ABSOLUTE-STATE CANCELLATION HYPOTHESIS SUPPORTED"
    else:
        verdict="TARGET STEP SOLVED WITH PHASE-RATE INITIALIZATION; REPRESENTATION-ONLY ROOT CAUSE NOT ISOLATED"
    timing={str(p.parent.relative_to(ROOT)):read(p) for p in ROOT.rglob("timing.json") if "tests" not in p.parts}
    summary={"moderate_accuracy":moderate["moderate_accuracy"],"target_tiny_step":target["target_tiny_step"],
        "cancellation_hypothesis":target["cancellation_hypothesis"],"scientific_verdict":verdict,
        "overall_verdict":"MODEL NOT YET VALIDATED","pilot_authorized":False,"full_series_authorized":False,
        "full_CHNS_transfer":"not_qualified","scientific_wall_s":sum(v["wall_s"] for v in timing.values()),
        "wall_budget_s":POLICY["wall_budget_s"],"timing":timing,
        "not_run":["50-step pilot","1068-step level0","half/quarter temporal series","full CHNS rate solver",
            "coupling probe","moving-contact benchmarks","tank/film/gravity/water-air"],
        "policy_sha256":json_hash(POLICY)}
    write(ROOT/"summary.json",summary)
    with (ROOT/"summary.csv").open("w") as stream:
        writer=csv.writer(stream); writer.writerow(["run","accepted","reason","iterations","initial_residual","final_residual","physical_checks"])
        for key,value in {**moderate.get("runs",{}),**target.get("runs",{})}.items():
            snes=value.get("SNES") or {}
            writer.writerow([key,value.get("accepted"),snes.get("reason"),snes.get("iterations"),
                snes.get("initial_residual"),snes.get("final_residual"),value.get("physical_checks")])
    print(json.dumps({k:v for k,v in summary.items() if k!="timing"}),flush=True)
    return summary


def all_stages():
    freeze(); moderate_accuracy()
    if read(ROOT/"moderate_accuracy/comparison.json")["moderate_accuracy"]=="qualified":
        if linear_tiny()["linear_comparison"]=="passed": target_tiny()
    return analyze()


STAGES={"freeze":freeze,"moderate-accuracy":moderate_accuracy,"linear-tiny":linear_tiny,
        "target-tiny":target_tiny,"analyze":analyze,"all":all_stages}
