"""STEP3A7 bounded stages. No new spectral analysis, optimizer or full CHNS."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import time
import xml.etree.ElementTree as ET
import numpy as np
from ..linearized_ch import array_hash, json_hash
from ..phase_rate_schedule import RateSchedule,cost_projection
from ..phase_rate_history import (read_json,atomic_json,file_hash,load_archive,Journal,
    validate_complete)
from ..nonlinear_accuracy import NonlinearControls
from ..provenance import source_hash
from ..phase_rate import VERSION
from .phase_rate_ch import HASHES,INITIAL_D
from . import spectral_ch as historical

ROOT=Path("validation_results/step3a7")
PLAN=Path("validation_results/step3a4/schedule_optimization/optimized_plan.json")
PRIOR_POLICY=Path("validation_results/step3a4/policy/policy.json")
SOURCE=Path(__file__).resolve().parents[1]
RUNNER_FILES=[Path("scripts/step3a7_benchmarks.py"),Path("scripts/step3a7_report.py"),Path("docker/step3/run.sh")]


def schedules():
    plan=read_json(PLAN)
    return [RateSchedule(plan,f,expected_hash=HASHES["optimized_plan_sha256"],expected_count=1068) for f in (1,2,4)]


def implementation():
    digest=hashlib.sha256()
    for path in RUNNER_FILES:
        digest.update(str(path).encode()); digest.update(path.read_bytes())
    return {"execution_base_HEAD":os.environ.get("STEP3_GIT_COMMIT","unavailable"),
        "multiphase_source_sha256":source_hash(),"runner_source_sha256":digest.hexdigest(),
        "docker_image":os.environ.get("STEP3_DOCKER_IMAGE_DIGEST","unavailable"),"phase_rate_version":VERSION}


def stack():
    import dolfinx, scipy, ufl
    from petsc4py import PETSc
    from mpi4py import MPI
    header=Path("/usr/local/petsc/linux-gnu-real64-32/include/dmumps_c.h").read_text()
    return {"DOLFINx":dolfinx.__version__,"PETSc":list(PETSc.Sys.getVersion()),
        "MUMPS":re.search(r'#define MUMPS_VERSION\s+"([^"]+)"',header).group(1),
        "scipy":scipy.__version__,"UFL":ufl.__version__,"scalar_precision":str(np.dtype(PETSc.ScalarType)),
        "precision_bits":np.finfo(PETSc.ScalarType).bits,"rank_count":MPI.COMM_WORLD.size}


def freeze():
    if (ROOT/"policy/frozen.json").exists():
        return source_guard()
    for filename in ("pure_preflight.xml","docker_preflight.xml"):
        node=ET.parse(ROOT/"tests"/filename).getroot()
        suites=list(node.iter("testsuite"))
        if not suites or any(int(s.get("failures",0))+int(s.get("errors",0)) for s in suites):
            raise ValueError("Passing pure and pinned tiny tests required BEFORE source freeze")
    prior=read_json(PRIOR_POLICY)
    if json_hash(prior["policy"])!=prior["policy_sha256"]:
        raise ValueError("Historical machine policy fingerprint mismatch")
    sparse=read_json("validation_results/step3a4/linear_sparse/optimized_plan/status.json")
    if sparse["qualification_status"]!="passed" or sparse["plan_sha256"]!=HASHES["optimized_plan_sha256"]:
        raise ValueError("Matching historical full-sparse PASS required")
    from ..equilibrium import equilibrium_config_fingerprint
    from ..config import ModelConfig
    prepared=read_json(historical.EQUILIBRIUM/"prepared.json")
    c=ModelConfig(**prepared["config"])
    policy={"inherited":prior["policy"],"inherited_policy_sha256":prior["policy_sha256"],
        "controls":NonlinearControls().record(),"hashes":HASHES,"version":"step3a7-v1",
        "first_guess":"semidiscrete","following_guess":"previous accepted phase_rate, no rescaling",
        "cost_statistic":"1.25*max(mean,p95) measured complete accepted intervals",
        "restart_field_relative":1e-10,"restart_D_relative":1e-10,"restart_energy_absolute":1e-12,
        "expected_wall_crossings":2,"physical_config_fingerprint":equilibrium_config_fingerprint(c),
        "target_mass":prepared["target_mass"],"initial_D_CH":INITIAL_D,"initial_certified_cells":8.328131075834218,
        "full_CHNS_authorized":False,"spectrum_optimizer_authorized":False}
    frozen={"implementation":implementation(),"stack":stack(),"policy":policy,"policy_sha256":json_hash(policy),
        "design_sha256":file_hash("STEP3A7_DESIGN.md"),
        "historical_input_provenance":read_json(historical.RESULTS/"linear_operator/status.json")["provenance"]}
    if frozen["stack"]["rank_count"]!=1 or frozen["implementation"]["execution_base_HEAD"]=="unavailable":
        raise ValueError("Pinned serial provenance required")
    (ROOT/"policy").mkdir(parents=True,exist_ok=True)
    with tarfile.open(ROOT/"policy/source.tar.gz","x:gz") as tar:
        for path in sorted(SOURCE.rglob("*.py")): tar.add(path,arcname="multiphase/"+str(path.relative_to(SOURCE)))
        for path in RUNNER_FILES: tar.add(path,arcname=str(path))
        tar.add("STEP3A7_DESIGN.md",arcname="STEP3A7_DESIGN.md")
    frozen["source_archive_sha256"]=file_hash(ROOT/"policy/source.tar.gz")
    for i,s in enumerate(schedules()): atomic_json(ROOT/"policy"/f"level{i}_schedule.json",s.descriptor)
    atomic_json(ROOT/"policy/frozen.json",frozen)
    return frozen


def source_guard():
    frozen=read_json(ROOT/"policy/frozen.json")
    if (frozen["implementation"]!=implementation() or frozen["stack"]!=stack() or
            frozen["design_sha256"]!=file_hash("STEP3A7_DESIGN.md") or
            frozen["policy_sha256"]!=json_hash(frozen["policy"]) or
            frozen["policy"]["inherited"]!=read_json(PRIOR_POLICY)["policy"]):
        raise ValueError("STOP scientific series: frozen source/policy/stack changed")
    for i,s in enumerate(schedules()):
        if read_json(ROOT/"policy"/f"level{i}_schedule.json")!=s.descriptor:
            raise ValueError("Frozen schedule changed")
    return frozen


def input_engine():
    from dolfinx import fem
    from ..phase_rate import IsolatedCHPhaseRate
    from ..equilibrium import equilibrium_config_fingerprint
    s,eq,psi,initial,provenance=historical.scientific_input()
    frozen=source_guard()
    if (array_hash(initial)!=HASHES["initial_phi_sha256"] or
            array_hash(psi)!=HASHES["perturbation_coefficient_sha256"] or
            provenance["equilibrium_fingerprint"]!=HASHES["equilibrium_fingerprint"] or
            equilibrium_config_fingerprint(s.config)!=frozen["policy"]["physical_config_fingerprint"]):
        raise ValueError("Historical physical input mismatch")
    field=fem.Function(eq.function_space); field.x.array[:]=initial; field.x.scatter_forward()
    engine=IsolatedCHPhaseRate(s,field)
    if array_hash(s.state.x.array[engine.physical_phi_map])!=HASHES["initial_phi_sha256"]:
        raise ValueError("Initialized phi hash mismatch")
    return engine,provenance


def run_identity(level, provenance):
    frozen=source_guard(); schedule=schedules()[level]
    return {"implementation":frozen["implementation"],"stack":frozen["stack"],
        "policy_sha256":frozen["policy_sha256"],"physical_config_fingerprint":frozen["policy"]["physical_config_fingerprint"],
        "mesh_fingerprint":provenance["mesh_fingerprint"],"input_hashes":HASHES,
        "schedule_sha256":schedule.sha256,"parent_plan_sha256":schedule.parent_sha256,
        "level":level,"rank_count":provenance["rank_count"]}


def sessions():
    return [read_json(p) for p in ROOT.glob("**/sessions/*/timing.json")]


def used_wall(level=None):
    return sum(r["wall_s"] for r in sessions() if level is None or r["level"]==level)


def execute(level, target, *, folder=None, resume=False, role="primary"):
    from ..phase_rate_trajectory import RateTrajectory,CostStop
    from ..performance import petsc_events
    from petsc4py import PETSc
    frozen=source_guard(); limits=frozen["policy"]["inherited"]
    folder=ROOT/f"level{level}" if folder is None else Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    session_root=folder/"sessions"; session_root.mkdir(exist_ok=True)
    session=session_root/f"{len(list(session_root.iterdir())):03d}"
    session.mkdir(exist_ok=False)
    atomic_json(session/"RUNNING.json",{"level":level,"target":target,"role":role,"identity":frozen["implementation"]})
    started=time.perf_counter(); runner=None; status=None; initialization=None
    PETSc.Log.begin()
    try:
        engine,provenance=input_engine(); identity=run_identity(level,provenance)
        runner=RateTrajectory(engine,schedules()[level],folder,identity,limits,frozen["policy"]["target_mass"],
            resume=resume,source_guard=source_guard)
        initialization=time.perf_counter()-started
        atomic_json(session/"provenance.json",{"identity":identity,"historical_input_provenance":frozen["historical_input_provenance"],
            "actual_input_provenance":provenance,"role":role})
        initial=runner.journal.rows[0]
        if (abs(initial["CH_dissipation"]/INITIAL_D-1)>1e-10 or
                abs(initial["cells_across_transition_certified_min"]/frozen["policy"]["initial_certified_cells"]-1)>1e-12):
            raise ValueError("Initial D/resolution reproduction mismatch; no PDE authorized")
        remaining=limits["isolated_series_wall_s"]-used_wall()-initialization
        if level==0: remaining=min(remaining,limits["isolated_level0_wall_s"]-used_wall(0)-initialization)
        if remaining<=0: raise CostStop("Frozen wall budget exhausted before trajectory")
        status=runner.run_until(target,wall_remaining_s=remaining)
        atomic_json(folder/"status.json",status)
        return status
    finally:
        if runner is not None:
            status=runner.status(); atomic_json(folder/"status.json",status); runner.close()
        atomic_json(session/"timing.json",{"level":level,"role":role,"wall_s":time.perf_counter()-started,
            "initialization_JIT_s":initialization,"target_steps":target,
            "accepted_steps":status["accepted_steps"] if status else None,"PETSc":petsc_events()})


def require_prefix(count):
    path=ROOT/"level0"
    if (path/"FAILED.json").exists(): raise ValueError("Scientific prefix failure; no downstream PDE")
    rows=Journal.read(path/"scalar_history.jsonl")
    if len(rows)<count+1 or not all(all(r["physical_gates"].values()) for r in rows):
        raise ValueError("Required accepted prefix missing")
    return rows


def pilot():
    try:
        status=execute(0,50)
    finally:
        if (ROOT/"level0/status.json").exists():
            status=read_json(ROOT/"level0/status.json")
            rows=Journal.read(ROOT/"level0/scalar_history.jsonl")
            atomic_json(ROOT/"pilot/status.json",{**status,"first_block_crossed":False,"restart_check":"not_run",
                "final_pilot_time":rows[-1]["time"],
                "max_Newton":max((r["solver"]["snes_iterations"] for r in rows[1:]),default=None),
                "max_weak_work_defect":max((abs(r["work"]["weak_work_defect"]) for r in rows[1:]),default=None),
                "level0_authorized":False,"dataset":"../level0 (reusable exact prefix)","temporal_qualification":False})
    return status


def restart_check():
    require_prefix(50); frozen=source_guard()
    src=ROOT/"level0"; dst=ROOT/"pilot/restart_fork"
    dst.mkdir(parents=True,exist_ok=False)
    arrays,meta=load_archive(src/"checkpoints/000025")
    shutil.copytree(src/"checkpoints/000025",dst/"checkpoints/000025")
    with (src/"scalar_history.jsonl").open("rb") as stream: prefix=stream.read(meta["journal_prefix"]["bytes"])
    with (dst/"scalar_history.jsonl").open("xb") as stream: stream.write(prefix)
    atomic_json(dst/"LATEST.json",{"checkpoint":"checkpoints/000025","checkpoint_sha256":meta["semantic_sha256"]})
    execute(0,28,folder=dst,resume=True,role="restart diagnostic fork; never replaces primary")
    from scipy.sparse import load_npz
    M=load_npz(historical.RESULTS/"linear_operator/M0.npz")
    def norm(v): return float(np.sqrt(max(0.,v @ (M@v))))
    comparisons=[]
    arows=Journal.read(src/"scalar_history.jsonl"); brows=Journal.read(dst/"scalar_history.jsonl")
    for i in (26,27,28):
        a,_=load_archive(src/"fields"/f"{i:06d}"); b,_=load_archive(dst/"fields"/f"{i:06d}")
        errors={key:norm(a[key]-b[key])/max(norm(a[key]),1e-30) for key in ("phi","mu","phase_rate")}
        errors.update(D=abs(arows[i]["cumulative_CH_dissipation"]-brows[i]["cumulative_CH_dissipation"])/max(abs(arows[i]["cumulative_CH_dissipation"]),1e-30),
            F=abs(arows[i]["E_total"]-brows[i]["E_total"]),budget=abs(arows[i]["energy_budget_defect"]-brows[i]["energy_budget_defect"]))
        passed=all(errors[k]<=frozen["policy"]["restart_field_relative"] for k in ("phi","mu","phase_rate"))
        passed &= errors["D"]<=frozen["policy"]["restart_D_relative"] and max(errors["F"],errors["budget"])<=frozen["policy"]["restart_energy_absolute"]
        comparisons.append({"step":i,"errors":errors,"passed":bool(passed)})
    result={"restart_check":"passed" if all(r["passed"] for r in comparisons) else "failed","comparisons":comparisons}
    atomic_json(ROOT/"pilot/restart_check.json",result)
    if result["restart_check"]!="passed": raise ValueError("Production restart comparison FAIL")
    return result


def forecast(next_level=0):
    frozen=source_guard(); all_s=schedules(); rows=[]; counts=[]
    for i,s in enumerate(all_s):
        path=ROOT/f"level{i}"/"timing_history.jsonl"
        if path.exists(): rows.extend(Journal.read(path))
        complete=ROOT/f"level{i}"/"COMPLETE.json"
        n=0
        if (ROOT/f"level{i}"/"scalar_history.jsonl").exists():
            n=Journal.read(ROOT/f"level{i}"/"scalar_history.jsonl")[-1]["step"]
        counts.append(0 if complete.exists() else s.nsteps-n)
    initialization=max(r["initialization_JIT_s"] for r in sessions() if r["initialization_JIT_s"] is not None)
    result=cost_projection(rows,counts,used_wall(),initialization,frozen["policy"]["inherited"])
    result["projected_level0_wall_s"]=used_wall(0)+result["remaining_wall_s"][0]
    result["level0_cost_ok"]=result["projected_level0_wall_s"]<=frozen["policy"]["inherited"]["isolated_level0_wall_s"]
    result["authorized"]=result["series_cost_ok"] and result["level0_cost_ok"]
    atomic_json(ROOT/"policy"/f"cost_before_level{next_level}.json",result)
    return result


def block_transition():
    require_prefix(50)
    if read_json(ROOT/"pilot/restart_check.json")["restart_check"]!="passed": raise ValueError("Restart gate required")
    schedule=schedules()[0]; status=execute(0,schedule.pilot_end,resume=True)
    rows=require_prefix(schedule.pilot_end); end=schedule.block_ends[0]
    counts=[r["solver"]["snes_iterations"] for r in rows[end-4:end+1]]
    post=rows[end+1]["solver"]
    good=post["snes_iterations"]<=max(5,2*max(counts))
    boundary={"first_block_end":end,"first_post_boundary":end+1,"before_iterations":counts,
        "post_boundary":post,"passed":good,"previous_rate_not_rescaled":post["guess"]["sha256"]==rows[end]["solver"]["retained_rate_sha256"]}
    boundary["passed"] &= boundary["previous_rate_not_rescaled"]
    atomic_json(ROOT/"block_transition/status.json",boundary)
    if not boundary["passed"]: raise ValueError("First dt-boundary gate failed")
    cost=forecast(0)
    result={**status,"first_block_crossed":True,"restart_check":"passed","physical_checks":"passed",
        "cost_projection":cost,"level0_authorized":cost["authorized"],"final_pilot_time":rows[-1]["time"],
        "max_Newton":max(r["solver"]["snes_iterations"] for r in rows[1:]),
        "max_weak_work_defect":max(abs(r["work"]["weak_work_defect"]) for r in rows[1:]),
        "dataset":"../level0 (continuous reusable prefix)","temporal_qualification":False}
    atomic_json(ROOT/"pilot/status.json",result)
    return result


def level(i):
    from ..phase_rate_trajectory import CostStop
    frozen=source_guard(); ss=schedules()
    if i==0:
        if not read_json(ROOT/"pilot/status.json")["level0_authorized"]: raise CostStop("Pilot did not authorize L0")
    else:
        for j in range(i):
            folder=ROOT/f"level{j}"; marker=read_json(folder/"COMPLETE.json")
            validate_complete(folder,marker["identity"],ss[j])
    cost=forecast(i)
    if not cost["authorized"]: raise CostStop("Frozen isolated series forecast gate failed")
    return execute(i,ss[i].nsteps,resume=(i==0))


def all_stages():
    freeze(); pilot(); restart_check(); block_transition()
    for i in range(3): level(i)
