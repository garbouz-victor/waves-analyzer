"""STEP3A9 production specialization of the STEP3A8 bounded stage harness.

Fresh outputs and conditional pressure-projection policy; solver core unchanged.
"""
import gc
import hashlib
import os
from pathlib import Path
import shutil
import tarfile
import time
import xml.etree.ElementTree as ET
import numpy as np
from ..linearized_ch import array_hash,json_hash
from ..phase_rate_history import (atomic_json,read_json,file_hash,archive,load_archive,
    Journal,validate_complete)
from ..full_rate_policy import CONTROLS,L0_SHA,cost_forecast as inherited_cost_forecast,session_result
from ..divergence_policy import production_policy,spatial_decision,bridge_divergence
from ..full_phase_rate import VERSION,CHNSPhaseRateBE
from ..full_rate_comparison import references
from ..phase_rate_trajectory import CostStop,ScientificFailure,rejected_metadata
from ..provenance import source_hash
from .phase_rate_series import schedules,stack
from .phase_rate_ch import HASHES
from . import spectral_ch as historical

ROOT=Path("validation_results/step3a9")
ISO=Path("validation_results/step3a7")
SOURCE=Path(__file__).resolve().parents[1]
RUNNERS=[Path("scripts/step3a9_transfer.py"),Path("scripts/step3a9_transfer_report.py"),Path("docker/step3/run.sh")]
POLICY=production_policy(read_json(ROOT/"divergence_audit/decision.json"))


def implementation():
    from .divergence_ch import implementation as audit_implementation
    core=audit_implementation()["solver_core_sha256"]
    h=hashlib.sha256()
    for p in RUNNERS: h.update(str(p).encode()); h.update(p.read_bytes())
    return {"execution_base_HEAD":os.environ.get("STEP3_GIT_COMMIT","unavailable"),
        "multiphase_source_sha256":source_hash(),"runner_sha256":h.hexdigest(),
        "original_weak_form_sha256":file_hash(SOURCE/"weak_form.py"),"full_phase_rate_version":VERSION,
        "docker_image":os.environ.get("STEP3_DOCKER_IMAGE_DIGEST","unavailable"),
        "solver_core_sha256":core,
        "diagnostics_sha256":file_hash(SOURCE/"divergence_observer.py"),
        "projection_sha256":file_hash(SOURCE/"divergence_audit.py")}


def audit_proof():
    """Verify the completed first freeze without pretending new files existed then."""
    f=read_json(ROOT/"policy/audit_frozen.json")
    if file_hash(ROOT/"policy/audit_source.tar.gz")!=f["source_archive_sha256"]:
        raise ScientificFailure("Audit archive changed")
    with tarfile.open(ROOT/"policy/audit_source.tar.gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():continue
            p=SOURCE/member.name[len("multiphase/"):] if member.name.startswith("multiphase/") else Path(member.name)
            if p.read_bytes()!=tar.extractfile(member).read():
                raise ScientificFailure("Frozen audit implementation changed: "+str(p))
    rows=[read_json(ROOT/"divergence_audit"/f"mesh_M{i}/summary.json") for i in range(3)]
    d=spatial_decision(rows);saved=read_json(ROOT/"divergence_audit/decision.json")
    if not d["passed"] or any(saved[k]!=v for k,v in d.items()):raise ScientificFailure("Audit decision not reproducible")
    accepted=[]
    for i in range(3):
        folder=ROOT/"divergence_audit"/f"mesh_M{i}/absolute_step"
        a,meta=load_archive(folder/"accepted");r=read_json(folder/"status.json")
        if not r["accepted"] or r["candidate_audit"]["physical"]["strong_divergence_L2"]!=rows[i]["strong_divergence_L2"]:
            raise ScientificFailure("Spatial measurement/candidate identity mismatch")
        accepted.append(meta["semantic_sha256"])
    sessions=[read_json(p) for p in ROOT.glob("divergence_audit/**/session_status.json")]
    if any(r["execution_status"]!="complete" for r in sessions):raise ScientificFailure("Incomplete audit session")
    return {"decision":saved,"accepted_archive_sha256":accepted,"audit_source_archive_sha256":f["source_archive_sha256"],
        "audit_wall_s":sum(r["wall_s"] for r in sessions),"stack":f["stack"]}


def audit_references():
    all_s=schedules()
    if all_s[0].sha256!=L0_SHA or all_s[0].nsteps!=1068: raise ValueError("L0 identity/count mismatch")
    refs=references(ISO,all_s)
    for ref in refs:
        if ref["identity"]["input_hashes"]!=HASHES: raise ValueError("Historical input hash mismatch")
    from scipy.sparse import load_npz
    M=load_npz(historical.RESULTS/"linear_operator/M0.npz")
    norm=lambda a:float(np.sqrt(max(0.,a@(M@a))))
    gaps=[]
    for n in all_s[0].common_parent_indices:
        a,_=load_archive(ISO/"level0/fields"/f"{n:06d}")
        b,_=load_archive(ISO/"level2/fields"/f"{4*n:06d}")
        gaps.append(dict(parent_step=n,phi_02_L2=norm(a["phi"]-b["phi"])))
    final0,final2=refs[0]["rows"][-1],refs[2]["rows"][-1]
    result={"passed":True,"COMPLETE":[r["marker"] for r in refs],"schedule_sha256":L0_SHA,
        "max_phi_02_L2":max(r["phi_02_L2"] for r in gaps),"phi_common_gaps":gaps,
        "integrated_D_02":abs(final0["cumulative_CH_dissipation"]-final2["cumulative_CH_dissipation"]),
        "final_F_02":abs(final0["E_total"]-final2["E_total"]),"history_and_field_hashes_verified":True}
    atomic_json(ROOT/"policy/isolated_references.json",result)
    return result


def tiny_bridge():
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..full_rate_bridge import audit
    c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(nx=4,nz=4,epsilon=1.,
        quadrature_degree=12,theta_equilibrium_deg=60.,g=0.,a_x=0.)
    from ..provenance import run_provenance
    s=CHNSSolver(c)
    result=audit(s); result.update(implementation=implementation(),stack=stack(),
        tiny_input_provenance=run_provenance(s),physical_config=c.to_dict() if hasattr(c,"to_dict") else vars(c),
        policy_sha256=json_hash(POLICY),scope="tiny arbitrary-field preflight, not production inputs")
    atomic_json(ROOT/"algebraic_bridge/tiny_fem/status.json",result)
    if not result["passed"]: raise ScientificFailure("Tiny full algebraic bridge failed")
    return result


def freeze():
    if (ROOT/"policy/production_frozen.json").exists(): return source_guard()
    for filename in ("production_pure_preflight.xml","production_docker_preflight.xml"):
        suites=list(ET.parse(ROOT/"tests"/filename).getroot().iter("testsuite"))
        if not suites or sum(int(s.get("tests",0)) for s in suites)==0 or any(int(s.get("failures",0))+int(s.get("errors",0)) for s in suites):
            raise ValueError("Passing host and pinned tiny tests required before freeze")
    proof=audit_proof()
    if proof["stack"]!=stack():raise ScientificFailure("Stack changed after spatial audit")
    bridge=read_json(ROOT/"tests/tiny_fem.json")
    if not bridge["passed"]: raise ValueError("Inherited unchanged-core preflight failed")
    refs=audit_references()
    inherited=read_json("validation_results/step3a4/policy/policy.json")
    if json_hash(inherited["policy"])!=inherited["policy_sha256"]: raise ValueError("Inherited policy corrupted")
    prepared=read_json(historical.EQUILIBRIUM/"prepared.json")
    prior=read_json(ISO/"policy/frozen.json")
    frozen={"implementation":implementation(),"stack":stack(),"policy":POLICY,"policy_sha256":json_hash(POLICY),
        "design_sha256":file_hash("STEP3A9_DESIGN.md"),"inherited_policy":inherited,
        "historical_input_provenance":read_json(historical.RESULTS/"linear_operator/status.json")["provenance"],
        "physical_config_fingerprint":prior["policy"]["physical_config_fingerprint"],
        "target_mass":prepared["target_mass"],"input_hashes":HASHES,"L0_schedule":schedules()[0].descriptor,
        "isolated_reference_markers":refs["COMPLETE"],"tiny_bridge_sha256":file_hash(ROOT/"tests/tiny_fem.json"),
        "audit_proof":proof,"audit_proof_sha256":json_hash(proof)}
    if frozen["stack"]["rank_count"]!=1 or frozen["implementation"]["execution_base_HEAD"]=="unavailable":
        raise ValueError("Pinned serial execution provenance required")
    with tarfile.open(ROOT/"policy/production_source.tar.gz","x:gz") as tar:
        for p in sorted(SOURCE.rglob("*.py")): tar.add(p,arcname="multiphase/"+str(p.relative_to(SOURCE)))
        for p in RUNNERS: tar.add(p,arcname=str(p))
        tar.add("STEP3A9_DESIGN.md",arcname="STEP3A9_DESIGN.md")
    frozen["source_archive_sha256"]=file_hash(ROOT/"policy/production_source.tar.gz")
    atomic_json(ROOT/"policy/production_policy.json",{"policy":POLICY,"policy_sha256":json_hash(POLICY),
        "audit_decision_sha256":file_hash(ROOT/"divergence_audit/decision.json")})
    atomic_json(ROOT/"policy/production_frozen.json",frozen)
    return frozen


def source_guard():
    f=read_json(ROOT/"policy/production_frozen.json")
    if (f["implementation"]!=implementation() or f["stack"]!=stack() or f["policy"]!=POLICY or
        f["policy_sha256"]!=json_hash(POLICY) or f["design_sha256"]!=file_hash("STEP3A9_DESIGN.md") or
        f["L0_schedule"]!=schedules()[0].descriptor or
        f["inherited_policy"]!=read_json("validation_results/step3a4/policy/policy.json") or
        f["source_archive_sha256"]!=file_hash(ROOT/"policy/production_source.tar.gz") or
        f["audit_proof_sha256"]!=json_hash(audit_proof())):
        raise ValueError("STOP: source/policy/stack/schedule changed after freeze")
    return f


def identity(provenance):
    f=source_guard()
    return {"implementation":f["implementation"],"source_archive_sha256":f["source_archive_sha256"],
        "stack":f["stack"],"policy_sha256":f["policy_sha256"],"input_hashes":HASHES,
        "physical_config_fingerprint":f["physical_config_fingerprint"],"mesh_fingerprint":provenance["mesh_fingerprint"],
        "schedule_sha256":L0_SHA,"level":"full_CHNS_L0","rank_count":f["stack"]["rank_count"]}


def input_engine(formulation="rate"):
    from dolfinx import fem
    from ..equilibrium import equilibrium_config_fingerprint
    from ..full_rate_experiment import AbsoluteBEReference
    frozen=source_guard(); s,eq,psi,initial,provenance=historical.scientific_input()
    if (array_hash(initial)!=HASHES["initial_phi_sha256"] or array_hash(psi)!=HASHES["perturbation_coefficient_sha256"] or
        provenance["equilibrium_fingerprint"]!=HASHES["equilibrium_fingerprint"] or
        equilibrium_config_fingerprint(s.config)!=frozen["physical_config_fingerprint"]):
        raise ValueError("Physical inputs mismatch; no PDE run")
    phi=fem.Function(eq.function_space); phi.x.array[:]=initial; phi.x.scatter_forward(); s.initialize(phi)
    e=CHNSPhaseRateBE(s) if formulation=="rate" else AbsoluteBEReference(s)
    if array_hash(s.state.x.array[e.physical_phi_map])!=HASHES["initial_phi_sha256"]:
        raise ValueError("Physical initial phi altered")
    return e,provenance


def sessions(): return [read_json(p) for p in ROOT.glob("**/sessions/*/session_status.json")]


def used(kind="total"):
    values=sessions()
    if any(r["status"]=="running" for r in values):
        raise ValueError("Unclosed scientific session: inspect crash before continuing")
    audit=read_json(ROOT/"policy/production_frozen.json")["audit_proof"]["audit_wall_s"] if kind=="total" else 0.
    return audit+sum(r["wall_s"] for r in values if kind=="total" or r.get(kind,False))


def cost_forecast(timings,remaining,used_full,used_total,initialization=0.):
    return inherited_cost_forecast(timings,remaining,used_full,used_total,initialization,policy=POLICY)


def start_session(folder,*,preliminary,full_path):
    source_guard(); total=used(); prelim=used("preliminary"); full=used("full_path")
    if total>=POLICY["total_wall_s"] or preliminary and prelim>=POLICY["preliminary_wall_s"] or full_path and full>=POLICY["full_L0_wall_s"]:
        raise CostStop("Frozen scientific wall budget exhausted")
    parent=Path(folder)/"sessions"; parent.mkdir(parents=True,exist_ok=True)
    path=parent/f"{len(list(parent.iterdir())):03d}"; path.mkdir(exist_ok=False)
    record={"status":"running","preliminary":preliminary,"full_path":full_path,"implementation":implementation()}
    atomic_json(path/"session_status.json",record)
    remaining=POLICY["total_wall_s"]-total
    if preliminary: remaining=min(remaining,POLICY["preliminary_wall_s"]-prelim)
    if full_path: remaining=min(remaining,POLICY["full_L0_wall_s"]-full)
    return path,record,remaining


def finish_session(path,record,started,error,**extra):
    atomic_json(path/"session_status.json",{**record,**extra,"status":session_result(error),
        "error":str(error) if error else None,"wall_s":time.perf_counter()-started})


def initial_check(row):
    if abs(row["CH_dissipation"]/POLICY["initial_D"]-1)>1e-10 or abs(row["cells_across_transition_certified_min"]/POLICY["initial_cells"]-1)>1e-12:
        raise ScientificFailure("Initial D_CH/certified resolution mismatch")
    eq=read_json(historical.EQUILIBRIUM/"prepared.json")
    # STEP3A7 supplies a directly measured initial excess; never infer a new one from a new solve.
    iso_initial=Journal.read(ISO/"level0/scalar_history.jsonl")[0]
    if abs(row["E_total"]-iso_initial["E_total"])>1e-13: raise ScientificFailure("Initial energy mismatch")
    if abs((row["F_CH"]-eq["free_energy"])-POLICY["initial_excess"])>1e-14:
        raise ScientificFailure("Initial excess energy mismatch")


def one_step(name,formulation,dt,level):
    from ..full_rate_experiment import solve_one
    from ..divergence_observer import DivergenceOneStepAudit as OneStepAudit
    from ..performance import petsc_events
    from petsc4py import PETSc
    folder=ROOT/name
    if folder.exists(): raise ValueError("Refusing repeated scientific one-step attempt: "+name)
    session,record,remaining=start_session(folder,preliminary=True,full_path=False)
    started=time.perf_counter(); error=None; initialization=None; PETSc.Log.begin()
    try:
        e,p=input_engine(formulation); ident=identity(p); f=source_guard()
        audit=OneStepAudit(e,f["inherited_policy"]["policy"],POLICY,f["target_mass"],production=True)
        initial_check(audit.initial); initialization=time.perf_counter()-started
        if initialization>=remaining: raise CostStop("Preliminary budget before solve")
        result=solve_one(e,dt,CONTROLS[level],audit,"step3a9_"+name.replace("/","_")+"_")
        result.update(identity=ident,formulation=formulation,accuracy_level=level,
            initial_guess="semidiscrete rate" if formulation=="rate" else "previous physical state",
            historical_input_provenance=f["historical_input_provenance"],actual_input_provenance=p)
        if result["accepted"]:
            archive_started=time.perf_counter()
            archive(folder/"accepted",e.last_snapshot,{"identity":ident,"accepted":True,"dt":dt,
                "formulation":formulation,"actual_controls":result["controls_after"]})
            result["archival_IO_s"]=time.perf_counter()-archive_started
        else:
            archive(folder/"rejected_attempt",{"last_residual_Function":e.state.x.array,
                "PETSc_iterate":e.problem.x.array,"last_accepted_physical_state":e.solver.state.x.array},
                {"identity":ident,"accepted":False,"execution_status":"failed","dt":dt})
        atomic_json(folder/"status.json",rejected_metadata(result))
        if not result["accepted"]: raise ScientificFailure(result["error"] or "Primary one-step not accepted")
        if time.perf_counter()-started>remaining: raise CostStop("Preliminary wall budget exceeded; no downstream work")
        return result,e.metric,audit.full
    except Exception as exc: error=exc; raise
    finally:
        finish_session(session,record,started,error,initialization_JIT_s=initialization,PETSc=petsc_events())


def moderate():
    source_guard()
    from ..full_rate_experiment import compare_states
    comparisons=[]
    for level in POLICY["moderate_levels"]:
        left=f"moderate/{level}_absolute"; right=f"moderate/{level}_rate"
        ar,_,_=one_step(left,"absolute",POLICY["moderate_dt"],level)
        gc.collect()
        br,metric,observer=one_step(right,"rate",POLICY["moderate_dt"],level)
        a,_=load_archive(ROOT/left/"accepted"); b,_=load_archive(ROOT/right/"accepted")
        c=compare_states(a,b,metric,observer,ar["candidate_audit"]["physical"],br["candidate_audit"]["physical"],POLICY)
        strong=bridge_divergence(ar["candidate_audit"]["physical"]["strong_divergence_L2"],
            br["candidate_audit"]["physical"]["strong_divergence_L2"])
        c["strong_divergence_same_root"]=strong
        c["passed"]=c["passed"] and strong["passed"]
        comparisons.append({"level":level,**c})
        atomic_json(ROOT/"moderate/comparison.json",{"comparisons":comparisons,
            "passed":len(comparisons)==len(POLICY["moderate_levels"]) and all(x["passed"] for x in comparisons),"temporal_qualification":False})
        print("moderate",level,c,flush=True); gc.collect()
    changes={}
    for formulation in ("absolute","rate"):
        a=read_json(ROOT/f"moderate/P1_{formulation}/status.json")["candidate_audit"]["physical"]["strong_divergence_L2"]
        b=read_json(ROOT/f"moderate/P2_{formulation}/status.json")["candidate_audit"]["physical"]["strong_divergence_L2"]
        changes[formulation]={"P1":a,"P2":b,"relative_change":abs(b-a)/max(a,b,1e-30),"diagnostic_only":True}
    atomic_json(ROOT/"moderate/Newton_stability.json",changes)
    if not all(x["passed"] for x in comparisons): raise ScientificFailure("Mandatory P1/P2 full formulation bridge unqualified")
    return comparisons


def target_tiny():
    if not read_json(ROOT/"moderate/comparison.json")["passed"]: raise ScientificFailure("Moderate prerequisite failed")
    r,metric,_=one_step("target_tiny/rate","rate",POLICY["tiny_dt"],"P0")
    from ..nonlinear_accuracy import read_snapshot
    iso,meta=read_snapshot("validation_results/step3a6/target_tiny/rate_semidiscrete")
    a,_=load_archive(ROOT/"target_tiny/rate/accepted")
    if float(iso["dt"])!=POLICY["tiny_dt"] or array_hash(iso["phi_old"])!=HASHES["initial_phi_sha256"]:
        raise ScientificFailure("Historical isolated tiny comparator mismatch")
    c={"accepted":r["accepted"],"dt":POLICY["tiny_dt"],"phi_L2_difference":metric.norm(a["phi_new"]-iso["phi_new"]),
        "mu_L2_difference":metric.norm(a["mu_new"]-iso["mu_new"]),
        "phi_relative_to_isolated":metric.norm(a["phi_new"]-iso["phi_new"])/metric.norm(iso["phi_new"]),
        "mu_relative_to_isolated":metric.norm(a["mu_new"]-iso["mu_new"])/metric.norm(iso["mu_new"]),
        "full_physical":r["candidate_audit"]["physical"],"full_work":r["candidate_audit"]["work"],
        "isolated_snapshot_sha256":meta["snapshot_sha256"],"coupling_qualification":False,
        "optional_absolute_control":"not_run; not elected in frozen design"}
    atomic_json(ROOT/"target_tiny/comparison.json",c); return c


def execute(target,*,folder=None,resume=False,preliminary=True,full_path=True):
    from ..divergence_trajectory import DivergenceTrajectory as FullRateTrajectory
    from ..performance import petsc_events
    from petsc4py import PETSc
    folder=ROOT/"full_coupling_L0" if folder is None else Path(folder)
    session,record,remaining=start_session(folder,preliminary=preliminary,full_path=full_path)
    started=time.perf_counter(); error=None; runner=None; initialization=None; PETSc.Log.begin()
    try:
        e,p=input_engine(); f=source_guard()
        runner=FullRateTrajectory(e,schedules()[0],folder,identity(p),f["inherited_policy"]["policy"],
            f["target_mass"],full_policy=POLICY,resume=resume,source_guard=source_guard)
        initial_check(runner.journal.rows[0]); initialization=time.perf_counter()-started
        if initialization>=remaining: raise CostStop("No remaining scientific wall budget")
        result=runner.run_until(target,wall_remaining_s=remaining-initialization)
        atomic_json(folder/"status.json",result); return result
    except Exception as exc: error=exc; raise
    finally:
        if runner is not None:
            atomic_json(folder/"status.json",runner.status()); runner.close()
        finish_session(session,record,started,error,initialization_JIT_s=initialization,target_steps=target,PETSc=petsc_events())


def require_prefix(n):
    folder=ROOT/"full_coupling_L0"
    if (folder/"FAILED.json").exists(): raise ScientificFailure("Failed primary forbids downstream PDE")
    rows=Journal.read(folder/"scalar_history.jsonl")
    if len(rows)<n+1 or not all(all(r["physical_gates"].values()) and all(r["full_physical_checks"].values()) for r in rows):
        raise ScientificFailure("Required fully checked continuous prefix missing")
    return rows


def pilot():
    if not read_json(ROOT/"target_tiny/comparison.json")["accepted"]: raise ScientificFailure("Target tiny prerequisite")
    try: return execute(50)
    finally:
        path=ROOT/"full_coupling_L0/status.json"
        if path.exists(): atomic_json(ROOT/"pilot/status.json",{**read_json(path),"first_block_crossed":False,
            "restart_check":"not_run","level0_authorized":False,"dataset":"../full_coupling_L0 (continuous reusable prefix)"})


def restart_check():
    require_prefix(50); source_guard()
    src=ROOT/"full_coupling_L0"; dst=ROOT/"pilot/restart_fork"; dst.mkdir(parents=True,exist_ok=False)
    a,meta=load_archive(src/"checkpoints/000025")
    shutil.copytree(src/"checkpoints/000025",dst/"checkpoints/000025")
    with (src/"scalar_history.jsonl").open("rb") as stream: prefix=stream.read(meta["journal_prefix"]["bytes"])
    with (dst/"scalar_history.jsonl").open("xb") as stream: stream.write(prefix)
    atomic_json(dst/"LATEST.json",{"checkpoint":"checkpoints/000025","checkpoint_sha256":meta["semantic_sha256"]})
    import subprocess,sys
    # Explicit new process, same pinned stack and source. No modified guess.
    subprocess.run([sys.executable,"scripts/step3a9_transfer.py","--stage","restart-fork"],check=True)
    # Reconstruct FE metrics only, not a PDE solve. Charge this read-only audit to preliminary wall.
    session,record,_=start_session(ROOT/"pilot/restart_audit",preliminary=True,full_path=False)
    started=time.perf_counter(); error=None
    try:
        from ..full_rate_diagnostics import FullObserver
        e,_=input_engine(); obs=FullObserver(e.solver,e.metric); result=[]
        left=Journal.read(src/"scalar_history.jsonl"); right=Journal.read(dst/"scalar_history.jsonl")
        for i in (26,27,28):
            x,_=load_archive(src/"fields"/f"{i:06d}"); y,_=load_archive(dst/"fields"/f"{i:06d}")
            errs={k:e.metric.norm(x[k]-y[k])/max(e.metric.norm(x[k]),1e-30) for k in ("phi","mu","phase_rate")}
            errs.update(u_scaled=obs.velocity_metric.norm(x["u"]-y["u"])/(obs.U_ref*np.sqrt(e.metric.area)),
                pi_scaled=obs.pressure_metric.norm(x["pi"]-y["pi"])/(obs.P_ref*np.sqrt(e.metric.area)))
            checks={k:bool(v<=(POLICY["restart_u_pi_scaled"] if k.endswith("scaled") else POLICY["restart_field_relative"])) for k,v in errs.items()}
            for k in ("cumulative_CH_dissipation","cumulative_viscous_dissipation","cumulative_slip_dissipation","E_total","energy_budget_defect"):
                gap=abs(left[i][k]-right[i][k]); relative=k=="cumulative_CH_dissipation"
                errs[k]=gap/max(abs(left[i][k]),1e-30) if relative else gap
                checks[k]=errs[k]<=(POLICY["restart_D_relative"] if relative else POLICY["restart_energy_absolute"])
            for k in ("strong_divergence_L2","projected_divergence_L2","orthogonal_divergence_L2","weak_continuity_Riesz"):
                gap=abs(left[i][k]-right[i][k])
                limit=POLICY["restart_field_relative"]*max(abs(left[i][k]),1e-30) if k in ("strong_divergence_L2","orthogonal_divergence_L2") else POLICY["weak_continuity"]
                errs[k]=gap;checks[k]=gap<=limit
            result.append({"step":i,"errors":errs,"checks":checks,"passed":all(checks.values())})
        out={"restart_check":"passed" if all(r["passed"] for r in result) else "failed","comparisons":result}
        atomic_json(ROOT/"pilot/restart_check.json",out)
        if out["restart_check"]!="passed": raise ScientificFailure("Production full restart mismatch")
        return out
    except Exception as exc: error=exc; raise
    finally: finish_session(session,record,started,error)


def block_transition():
    require_prefix(50)
    if read_json(ROOT/"pilot/restart_check.json")["restart_check"]!="passed": raise ScientificFailure("Restart prerequisite")
    schedule=schedules()[0]; execute(schedule.pilot_end,resume=True)
    rows=require_prefix(schedule.pilot_end); end=schedule.block_ends[0]
    before=[r["solver"]["snes_iterations"] for r in rows[end-4:end+1]]
    post=rows[end+1]["solver"]
    passed=post["reason"]>0 and post["snes_iterations"]<=max(5,2*max(before)) and post["guess"]["sha256"]==rows[end]["solver"]["retained_rate_sha256"]
    boundary={"first_block_end":end,"first_post_boundary":end+1,"before_iterations":before,"post":post,"passed":passed}
    atomic_json(ROOT/"pilot/block_transition.json",boundary)
    if not passed: raise ScientificFailure("First full dt-block transition unqualified")
    times=Journal.read(ROOT/"full_coupling_L0/timing_history.jsonl")
    cost=cost_forecast(times,schedule.nsteps-rows[-1]["step"],used("full_path"),used(),
        max(r.get("initialization_JIT_s",0.) or 0. for r in sessions()))
    atomic_json(ROOT/"policy/cost_before_full_L0.json",cost)
    atomic_json(ROOT/"pilot/status.json",{**read_json(ROOT/"full_coupling_L0/status.json"),
        "first_block_crossed":True,"restart_check":"passed","level0_authorized":cost["authorized"],"cost_projection":cost,
        "final_pilot_time":rows[-1]["time"],"max_Newton":max(r["solver"]["snes_iterations"] for r in rows[1:]),
        "max_weak_work":max(abs(r["work"]["weak_work_defect"]) for r in rows[1:])})
    return cost


def full_L0():
    if not read_json(ROOT/"pilot/status.json")["level0_authorized"]: raise CostStop("Pilot/cost gate does not authorize full horizon")
    require_prefix(schedules()[0].pilot_end)
    return execute(schedules()[0].nsteps,resume=True,preliminary=False)


def all_stages():
    freeze(); moderate(); target_tiny(); pilot(); restart_check(); block_transition(); full_L0()
