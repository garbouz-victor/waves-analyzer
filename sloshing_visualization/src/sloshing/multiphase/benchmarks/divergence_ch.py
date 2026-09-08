"""STEP3A9 PART A. All output is new; no historical candidate is accepted."""
import contextlib
import hashlib
import os
from pathlib import Path
import signal
import tarfile
import time
import xml.etree.ElementTree as ET
import numpy as np
from ..divergence_policy import AUDIT,spatial_decision
from ..phase_rate_history import read_json,atomic_json,archive,load_archive,file_hash
from ..linearized_ch import array_hash,json_hash
from ..provenance import source_hash,run_provenance
from ..phase_rate_trajectory import ScientificFailure,rejected_metadata
from . import spectral_ch as historical
from .phase_rate_series import stack,schedules
from .phase_rate_ch import HASHES

ROOT=Path("validation_results/step3a9")
SOURCE=Path(__file__).resolve().parents[1]
CORE=("full_phase_rate.py","full_rate_bridge.py","weak_form.py","free_energy.py","material.py","boundary.py")
RUNNERS=(Path("scripts/step3a9_benchmarks.py"),Path("scripts/step3a9_report.py"),Path("docker/step3/run.sh"))


def implementation():
    core={p:file_hash(SOURCE/p) for p in CORE}
    with tarfile.open("validation_results/step3a8/policy/source.tar.gz") as old:
        archived={p:hashlib.sha256(old.extractfile("multiphase/"+p).read()).hexdigest() for p in CORE}
    if core!=archived: raise ScientificFailure("Solver core changed from STEP3A8; STOP")
    h=hashlib.sha256()
    for p in RUNNERS:h.update(str(p).encode());h.update(p.read_bytes())
    return {"execution_base_HEAD":os.environ.get("STEP3_GIT_COMMIT","unavailable"),
        "multiphase_source_sha256":source_hash(),"solver_core_sha256":core,
        "runner_sha256":h.hexdigest(),"docker_image":os.environ.get("STEP3_DOCKER_IMAGE_DIGEST","unavailable")}


def verify_references():
    from scipy.sparse import load_npz
    from ..full_rate_comparison import references
    ss=schedules(); refs=references(Path("validation_results/step3a7"),ss)
    if ss[0].sha256!="21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31" or ss[0].nsteps!=1068:
        raise ValueError("L0 schedule mismatch")
    if any(r["identity"]["input_hashes"]!=HASHES for r in refs):raise ValueError("Historical inputs mismatch")
    M=load_npz(historical.RESULTS/"linear_operator/M0.npz")
    gaps=[]
    for n in ss[0].common_parent_indices:
        a,_=load_archive(Path("validation_results/step3a7/level0/fields")/f"{n:06d}")
        b,_=load_archive(Path("validation_results/step3a7/level2/fields")/f"{4*n:06d}")
        d=a["phi"]-b["phi"];gaps.append(float(np.sqrt(d@(M@d))))
    a,b=refs[0]["rows"][-1],refs[2]["rows"][-1]
    return {"verified_COMPLETE":[r["marker"] for r in refs],"schedule_sha256":ss[0].sha256,
        "max_phi_02_L2":max(gaps),"integrated_D_02":abs(a["cumulative_CH_dissipation"]-b["cumulative_CH_dissipation"]),
        "final_F_02":abs(a["E_total"]-b["E_total"]),"all_history_and_field_hashes_verified":True}


def freeze():
    path=ROOT/"policy/audit_frozen.json"
    if path.exists():return guard()
    for name in ("pure_preflight.xml","docker_preflight.xml"):
        suites=list(ET.parse(ROOT/"tests"/name).getroot().iter("testsuite"))
        if not suites or not sum(int(s.get("tests",0)) for s in suites) or any(int(s.get("failures",0))+int(s.get("errors",0)) for s in suites):
            raise ValueError("Passing pure and pinned preflight required")
    small=read_json(ROOT/"tests/tiny_fem.json")
    if not small["passed"] or small["implementation"]!=implementation():
        raise ValueError("Tiny projection/kernel/bridge evidence must match source")
    refs=verify_references()
    eq=read_json(historical.EQUILIBRIUM/"prepared.json")
    if [eq["config"]["nx"],eq["config"]["nz"]]!=AUDIT["meshes"][0]:raise ValueError("Production mesh mismatch")
    if eq["target_mass"]!=AUDIT["target_mass"] or eq["fingerprint"]!=HASHES["equilibrium_fingerprint"]:
        raise ValueError("Mass/equilibrium mismatch")
    inherited=read_json("validation_results/step3a4/policy/policy.json")
    if json_hash(inherited["policy"])!=inherited["policy_sha256"]:raise ValueError("Inherited policy SHA mismatch")
    old=Path("validation_results/step3a8/moderate/P1_absolute")
    _,meta=load_archive(old/"rejected_attempt")
    f={"implementation":implementation(),"stack":stack(),"audit_policy":AUDIT,
        "decision_rule_sha256":json_hash(AUDIT),"mesh_hierarchy_sha256":json_hash(AUDIT["meshes"]),
        "design_sha256":file_hash("STEP3A9_DESIGN.md"),"input_hashes":HASHES,
        "inherited_policy":inherited,"isolated_references":refs,
        "historical_STEP3A8_STOP":{"status":read_json(old/"status.json"),"rejected_archive_metadata":meta,
            "status_sha256":file_hash(old/"status.json"),"still_rejected":True},
        "historical_prepared_metadata":eq}
    if f["implementation"]["execution_base_HEAD"]!="5f2d9498e1d0715f57b125ee49b4ddfd533467ae" or f["stack"]["rank_count"]!=1:
        raise ValueError("Expected base/serial pinned execution required")
    (ROOT/"policy").mkdir(parents=True,exist_ok=True)
    with tarfile.open(ROOT/"policy/audit_source.tar.gz","x:gz") as tar:
        for p in sorted(SOURCE.rglob("*.py")):tar.add(p,arcname="multiphase/"+str(p.relative_to(SOURCE)))
        for p in RUNNERS:tar.add(p,arcname=str(p))
        tar.add("STEP3A9_DESIGN.md",arcname="STEP3A9_DESIGN.md")
    f["source_archive_sha256"]=file_hash(ROOT/"policy/audit_source.tar.gz")
    atomic_json(path,f);return f


def preflight():
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..divergence_audit import kernel_demonstration
    from ..full_rate_bridge import audit
    c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(nx=4,nz=4,epsilon=1.,
        quadrature_degree=12,theta_equilibrium_deg=60.,g=0.,a_x=0.)
    s=CHNSSolver(c);k=kernel_demonstration(s);bridge=audit(s)
    result={"kernel":k,"bridge":bridge,"passed":k["passed"] and bridge["passed"],
        "implementation":implementation(),"stack":stack(),"scope":"tiny mathematical FE preflight only"}
    atomic_json(ROOT/"tests/tiny_fem.json",result)
    if not result["passed"]:raise ScientificFailure("Tiny preflight failed")
    print(result,flush=True)


def guard():
    f=read_json(ROOT/"policy/audit_frozen.json")
    if (f["implementation"]!=implementation() or f["stack"]!=stack() or f["audit_policy"]!=AUDIT or
        f["decision_rule_sha256"]!=json_hash(AUDIT) or f["design_sha256"]!=file_hash("STEP3A9_DESIGN.md") or
        f["source_archive_sha256"]!=file_hash(ROOT/"policy/audit_source.tar.gz")):
        raise ScientificFailure("Audit source/policy/stack changed; STOP")
    return f


@contextlib.contextmanager
def session(folder):
    frozen=guard();folder=Path(folder)
    used=sum(r["wall_s"] for r in (read_json(p) for p in ROOT.glob("divergence_audit/**/session_status.json")))
    remaining=AUDIT["wall_budget_s"]-used
    if remaining<=0:raise ScientificFailure("Audit wall budget exhausted")
    if folder.exists():raise ValueError("No repeated scientific attempt: "+str(folder))
    folder.mkdir(parents=True)
    record={"execution_status":"running","implementation":frozen["implementation"],"stack":frozen["stack"],
        "source_archive_sha256":frozen["source_archive_sha256"],"decision_rule_sha256":frozen["decision_rule_sha256"],
        "wall_s":0.,"budget_remaining_s":remaining}
    atomic_json(folder/"session_status.json",record)
    started=time.perf_counter();error=None
    def timeout(*_):raise ScientificFailure("Frozen audit wall budget exceeded")
    previous=signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,remaining)
    try:
        yield record
        guard()
        if time.perf_counter()-started>remaining:raise ScientificFailure("Audit wall budget exceeded")
    except BaseException as exc:error=exc;raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous)
        record.update(execution_status="failed" if error else "complete",wall_s=time.perf_counter()-started,
                      error=str(error) if error else None)
        atomic_json(folder/"session_status.json",record)


def prepare_mesh(index):
    from dolfinx import fem
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..equilibrium import solve_equilibrium,zero_mass_perturbations,EquilibriumFailure
    from ..initialization import sessile_drop
    from ..diagnostics import Diagnostics,interface_resolution
    from ..ch_validation_diagnostics import WallTrace
    folder=ROOT/"divergence_audit"/f"mesh_M{index}"/"initialization"
    with session(folder) as record:
        eq_meta=read_json(historical.EQUILIBRIUM/"prepared.json")
        if index==0:
            s,eq,psi,initial,prov=historical.scientific_input();prepared_phi=eq.x.array.copy()
            if array_hash(initial)!=HASHES["initial_phi_sha256"] or array_hash(psi)!=HASHES["perturbation_coefficient_sha256"]:
                raise ScientificFailure("M0 input hash mismatch")
        else:
            nx,nz=AUDIT["meshes"][index]
            c=ModelConfig(**eq_meta["config"]).changed(nx=nx,nz=nz,max_unknowns=AUDIT["max_unknowns_refined"])
            s=CHNSSolver(c);s.initialize(sessile_drop(c.epsilon,c.refinement_circle_radius,60.,c.z_min))
            started=time.perf_counter()
            try:
                prepared=solve_equilibrium(s,target_mass=AUDIT["target_mass"],mu_initial=c.sigma/(2*c.refinement_circle_radius),
                    history_callback=lambda r:print(f"M{index} equilibrium "+str(r),flush=True))
            except EquilibriumFailure as exc:
                archive(folder/"rejected_equilibrium",{"phi":exc.coefficients},
                    {"accepted":False,"diagnostics":exc.diagnostics})
                raise
            record["equilibrium_s"]=time.perf_counter()-started
            eq_meta=prepared.metadata;s.initialize_prepared_equilibrium(prepared)
            prepared_phi=prepared.phi_coefficients.copy();eq=s.state.sub(2).collapse()
            psi=zero_mass_perturbations(s,eq)["smooth_bulk"][0].x.array.copy()
            initial=prepared_phi+AUDIT["delta"]*psi;prov=run_provenance(s)
        equilibrium_row=Diagnostics(s).measure();eq_res=interface_resolution(s)
        mu=s.state.sub(3).collapse().x.array
        eq_checks={"stationarity":eq_meta["stationarity"]["normalized"]<=1e-8,
            "mass":eq_meta["mass_residual_relative_to_domain"]<=1e-10,
            "constant_mu":float(np.ptp(mu))<=1e-12,
            "CH_zero":equilibrium_row["CH_dissipation"]<=1e-16,
            "resolved":eq_res["cells_across_transition_certified_min"]>=8}
        if not all(eq_checks.values()):raise ScientificFailure("Initial equilibrium qualification failed: "+str(eq_checks))
        phi=fem.Function(eq.function_space);phi.x.array[:]=initial;phi.x.scatter_forward();s.initialize(phi)
        row=Diagnostics(s).measure();resolution=interface_resolution(s)
        if resolution["cells_across_transition_certified_min"]<8 or len(WallTrace(s).crossings())!=2:
            raise ScientificFailure("Initial resolution/topology failed")
        if index==0 and (abs(row["CH_dissipation"]/.11400313905309394-1)>1e-10 or
                abs(resolution["cells_across_transition_certified_min"]/8.328131075834218-1)>1e-12):
            raise ScientificFailure("M0 initial D/certification mismatch")
        meta={"accepted":True,"role":"t=0 initialization, not a transient endpoint", "provenance":prov,
            "physical_config":s.config.as_dict(),"prepared":eq_meta,"equilibrium_checks":eq_checks,
            "initial_metrics":row,"initial_resolution":resolution,"mesh_report":s.mesh_report,
            "initial_construction":"historical exact M0; same analytic continuum family with own discrete equilibrium on M1/M2"}
        archive(folder/"initial",{"phi_eq":prepared_phi,"psi":psi,"initial_phi":initial,"physical_state":s.state.x.array},meta)
        atomic_json(folder/"status.json",{"passed":True,"equilibrium_checks":eq_checks,"initial_D":row["CH_dissipation"],
            "certified_cells":resolution["cells_across_transition_certified_min"]})


def rejected_audit():
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..phase_rate import PhaseMetric
    from ..full_rate_diagnostics import FullObserver
    from ..divergence_audit import DivergenceAudit
    from ..diagnostics import Diagnostics
    from ..full_rate_policy import POLICY
    folder=ROOT/"divergence_audit/rejected_step3a8"
    with session(folder) as record:
        old=Path("validation_results/step3a8/moderate/P1_absolute/rejected_attempt")
        arrays,meta=load_archive(old)
        s=CHNSSolver(ModelConfig(**read_json(historical.EQUILIBRIUM/"prepared.json")["config"]))
        s.state.x.array[:]=arrays["last_residual_Function"];s.state.x.scatter_forward()
        row=Diagnostics(s).measure();metric=PhaseMetric(s,s.space.sub(2).collapse()[0])
        existing=FullObserver(s,metric).measure(row,POLICY)
        audit=DivergenceAudit(s);result=audit.measure(existing["weak_continuity_Riesz"])
        cells,statistics=audit.cells()
        result.update(cellwise=statistics,historical_existing_observer=existing,
            original_archive_sha256=meta["semantic_sha256"],accepted=False,
            role="read-only historical REJECTED candidate; no solve or publication",
            PETSc_Function_equal=bool(np.array_equal(arrays["last_residual_Function"],arrays["PETSc_iterate"])))
        archive(folder/"projection",{"projected_divergence":audit.projected.x.array,**cells},
                {"accepted":False,"role":"read-only diagnostic","historical_archive":meta["semantic_sha256"]})
        atomic_json(folder/"status.json",result)
        plot_cells(folder,cells,result)
        if not result["projection_agrees"] or not result["pythagorean_pass"]:raise ScientificFailure("Independent projection audit failed")


def mesh_step(index):
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..full_rate_experiment import AbsoluteBEReference,solve_one
    from ..divergence_observer import DivergenceOneStepAudit
    from ..full_rate_policy import POLICY,CONTROLS
    for i in range(3):
        if not read_json(ROOT/"divergence_audit"/f"mesh_M{i}/initialization/status.json")["passed"]:
            raise ScientificFailure("All three initial states must qualify before CHNS")
    folder=ROOT/"divergence_audit"/f"mesh_M{index}"/"absolute_step"
    with session(folder) as record:
        arrays,meta=load_archive(folder.parent/"initialization/initial")
        s=CHNSSolver(ModelConfig(**meta["physical_config"]))
        for f in (s.state,s.old,s.older):f.x.array[:]=arrays["physical_state"];f.x.scatter_forward()
        s._check_material();e=AbsoluteBEReference(s)
        frozen=guard();audit=DivergenceOneStepAudit(e,frozen["inherited_policy"]["policy"],POLICY,AUDIT["target_mass"])
        result=solve_one(e,AUDIT["dt"],CONTROLS["P1"],audit,f"step3a9_spatial_M{index}_")
        result.update(actual_input_provenance=run_provenance(s),initialization_metadata_sha256=meta["semantic_sha256"],
            decision_rule_sha256=frozen["decision_rule_sha256"],implementation=frozen["implementation"],stack=frozen["stack"])
        if result["accepted"]:
            archive(folder/"accepted",e.last_snapshot,{"accepted":True,"role":"spatial audit only", "dt":AUDIT["dt"],
                "initialization_metadata_sha256":meta["semantic_sha256"],"actual_controls":result["controls_after"]})
            row=result["candidate_audit"]["physical"]
            tri=s.mesh.geometry.x[s.mesh.geometry.dofmap,:2]
            h=float(np.linalg.norm(tri-np.roll(tri,1,axis=1),axis=2).max())
            summary={**row,"h":h,"mesh_shape":AUDIT["meshes"][index],"cells":len(tri),
                "velocity_DOFs":len(e.u_map),"pressure_DOFs":len(e.pi_map),
                "physical_pass":all(result["candidate_audit"]["physical_gates"].values()) and all(row["full_physical_checks"].values()),
                "SNES_reason":result["SNES_reason"],"KSP_success":all(r["KSP_reason"]>0 for r in result["trace"][1:])}
            atomic_json(folder.parent/"summary.json",summary)
        else:
            archive(folder/"rejected_attempt",{"last_residual_Function":e.state.x.array,"PETSc_iterate":e.problem.x.array,
                "last_accepted_physical_state":s.state.x.array},{"accepted":False,"execution_status":"failed"})
        atomic_json(folder/"status.json",rejected_metadata(result))
        if not result["accepted"]:raise ScientificFailure("Spatial primary failed: "+str(result["error"]))


def decision():
    guard();rows=[read_json(ROOT/"divergence_audit"/f"mesh_M{i}/summary.json") for i in range(3)]
    d=spatial_decision(rows)
    atomic_json(ROOT/"divergence_audit/spatial_convergence.json",{"rows":rows,**d})
    atomic_json(ROOT/"divergence_audit/decision.json",{**d,"decision_rule_sha256":json_hash(AUDIT)})
    import matplotlib;matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots();ax.loglog([r["h"] for r in rows],[r["strong_divergence_L2"] for r in rows],"o-",label=f"measured slope {d['p_div']:.5g}")
    ax.set(xlabel="max cell diameter [m]",ylabel="strong divergence L2");ax.legend();fig.tight_layout()
    fig.savefig(ROOT/"divergence_audit/divergence_spatial_convergence.pdf");plt.close(fig)
    return d


def plot_cells(folder,arrays,result):
    import matplotlib;matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from matplotlib.colors import LogNorm
    xy=arrays["triangles"].reshape(-1,2);tri=Triangulation(xy[:,0],xy[:,1],np.arange(len(xy)).reshape(-1,3))
    fig,ax=plt.subplots(figsize=(9,4.7));values=arrays["contributions"]
    artist=ax.tripcolor(tri,facecolors=values,norm=LogNorm(vmin=max(float(values.min()),1e-40),vmax=float(values.max())),rasterized=True)
    ax.tricontour(tri,arrays["phi_vertices"].ravel(),levels=[0.],colors="white",linewidths=.8)
    ax.axhline(0.,color="black");ax.set(xlabel="x [m]",ylabel="z [m]",aspect="equal",
        title="Historical rejected state: cell integral of (div u)^2; vertex contour is visual only")
    fig.colorbar(artist,ax=ax,label="cell squared-divergence contribution");fig.tight_layout()
    fig.savefig(folder/"strong_divergence_cells.pdf");plt.close(fig)
    fig,ax=plt.subplots();keys=["strong_divergence_L2","projected_divergence_L2","orthogonal_divergence_L2"]
    ax.bar(["strong","P_Q divergence","orthogonal"],[result[k] for k in keys]);ax.set_yscale("log")
    ax.set_ylabel("L2 norm");fig.tight_layout();fig.savefig(folder/"divergence_projection.pdf");plt.close(fig)
