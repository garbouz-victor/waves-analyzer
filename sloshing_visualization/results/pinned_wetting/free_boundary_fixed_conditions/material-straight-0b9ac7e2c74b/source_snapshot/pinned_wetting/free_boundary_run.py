"""Small PW2 native/checkpoint adapter, using the existing single job ledger."""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import shutil
import time
from types import SimpleNamespace
import h5py
import numpy as np
from . import free_boundary as fb
from .coating import Coating
from .mission import doctor, read_json, write_json, sha256, job, lock_busy

REL_RUNS = "sloshing_visualization/results/pinned_wetting/free_boundary_fixed_conditions"
REL_OUTPUT = "sloshing_visualization/output/pinned_wetting/free_boundary_fixed_conditions"


def identity(root, c):
    contract = read_json(root/"mission/pinned_wetting/solve_to_animation_v2/CONTRACT.json")
    files = ["pinned_wetting/free_boundary.py", "pinned_wetting/free_boundary_run.py",
             "pinned_wetting/free_boundary_remesh.py", "pinned_wetting/coating.py", "mesh.py", "config.py"]
    hashes = {p:sha256(root/"sloshing_visualization/src/sloshing"/p) for p in files}
    versions = doctor(root)["versions"]
    numerical = {"model_id":fb.MODEL_ID, "physical_contract":contract["physical_contract"],
        "controls":asdict(c), "source_hashes":hashes,
        "libraries":{k:versions[k] for k in ("numpy","scipy","scikit-fem","h5py")},
        "method":"full_material_P2_P1_nonlinear_midpoint", "synthetic":False}
    return numerical, hashlib.sha256(json.dumps(numerical,sort_keys=True).encode()).hexdigest()


def save_csr(group, matrix):
    matrix = matrix.copy().tocsr(); matrix.sum_duplicates(); matrix.sort_indices()
    for key in ("data", "indices", "indptr"): group.create_dataset(key,data=getattr(matrix,key))
    group.attrs["shape"] = matrix.shape


def save_state(h, step, t, mesh, v, p, coating, row, rezoning=None):
    group = h["states"].create_group(f"{step:08d}")
    for key, value in (("geometry",mesh.p),("velocity",v),("pressure_midpoint",p)):
        group.create_dataset(key,data=value,compression="gzip",compression_opts=1)
    group.attrs.update(step=step,time_s=t,accepted=True,coating=json.dumps(coating.to_dict()),
                       diagnostics=json.dumps(row))
    if rezoning is not None:
        transfer=group.create_group("rezoning_before_step")
        for key,value in rezoning.items():
            if isinstance(value,np.ndarray): transfer.create_dataset(key,data=value,compression="gzip",compression_opts=1)
            else: transfer.attrs[key]=value
    h.attrs["last_accepted_step"] = step
    h.attrs["last_accepted_time_s"] = t
    h.flush()


def run_case(root, c, role, state, heartbeat, stop_after=None):
    numerical, digest = identity(root,c)
    run_id = f"material-{role}-{digest[:12]}"
    run_dir = root/REL_RUNS/run_id
    created = not run_dir.exists()
    run_dir.mkdir(parents=True,exist_ok=True)
    if created:
        write_json(run_dir/"identity.json",numerical)
        write_json(run_dir/"resolved_case.json",{"model_id":fb.MODEL_ID,
            "physical_contract":numerical["physical_contract"],"numerics":asdict(c),
            "provenance":{"actual_execution_HEAD":doctor(root)["execution_HEAD"],"run_id":run_id}})
        for name in numerical["source_hashes"]:
            dest = run_dir/"source_snapshot"/name; dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(root/"sloshing_visualization/src/sloshing"/name,dest)
        shutil.copyfile(root/"mission/pinned_wetting/solve_to_animation_v2/MODEL_NOTES.md",run_dir/"model_notes.md")
    elif read_json(run_dir/"identity.json") != numerical:
        raise ValueError("PW2 incompatible numerical identity")
    head = doctor(root)["execution_HEAD"]
    provenance = read_json(run_dir/"provenance.json") if (run_dir/"provenance.json").exists() else {
        "execution_HEAD_at_creation":head,"execution_HEAD_history":[]}
    if not provenance["execution_HEAD_history"] or provenance["execution_HEAD_history"][-1]["HEAD"]!=head:
        provenance["execution_HEAD_history"].append({"HEAD":head,"time":time.time(),"action":"create_or_resume"})
    write_json(run_dir/"provenance.json",provenance)
    state.update(current_run_id=run_id,current_run_scope="PW2_FREE_BOUNDARY_FIXED_CONDITIONS",
        current_model_id=fb.MODEL_ID,status="RUNNING",last_error=None)
    if run_id not in state.setdefault("run_ids",[]): state["run_ids"].append(run_id)
    state.setdefault("free_boundary_runs",{})[role] = {"run_id":run_id,"path":str(run_dir),"controls":asdict(c)}
    native = run_dir/"native.h5"
    mesh = fb.initial_mesh(c)
    op0 = fb.Operators(mesh,c)
    ids = fb.interface_nodes(mesh)
    _, orientation = fb.quadratic_minima(mesh)
    E0 = 9.81*np.tan(np.deg2rad(2.))**2/3.
    rejected = 0
    with h5py.File(native,"a") as h:
        if "states" not in h:
            h.attrs.update(model_id=fb.MODEL_ID,numerical_identity=json.dumps(numerical),
                numerical_digest=digest,run_id=run_id,synthetic=False,equations="full_material_NS",
                pressure_variable="pi=p+gz" if c.hydrostatic_split else "gauge_p")
            topo = h.create_group("topology")
            topo.create_dataset("triangles",data=mesh.t)
            topo.create_dataset("initial_geometry",data=mesh.p)
            topo.create_dataset("interface_nodes",data=ids)
            topo.create_dataset("orientation",data=orientation)
            for side, facets in mesh.boundaries.items(): topo.create_dataset(side+"_facets",data=facets)
            topo.create_dataset("fixed_velocity_dofs",data=op0.fixed)
            topo.create_dataset("component_dofs",data=np.array(op0.components))
            save_csr(h.create_group("initial_left_wall_operator"),op0.Kl)
            h.create_group("states"); h.create_group("rejected_trials")
            acceleration, p = op0.startup()
            h.create_dataset("initial_acceleration",data=acceleration)
            v = np.zeros(op0.vbasis.N)
            coating = Coating.initial(run_id,-10.,mesh.p[1,ids[0]],mesh.p[1,ids[-1]])
            row = op0.measure(v)
            row.update(cumulative_bulk_loss=0.,cumulative_left_wall_loss=0.,
                cumulative_geometry_kinetic_work=0.,cumulative_geometry_potential_work=0.,
                cumulative_rezone_work=0.,last_rezone_time_s=0.,
                continuous_energy_error=(row["energy"]-E0)/E0,dt_s=0.)
            save_state(h,0,0.,mesh,v,p,coating,row)
        # Only explicitly committed native states can be resumed.
        n = int(h.attrs["last_accepted_step"])
        prev = h["states"][f"{n:08d}"]
        mesh = replace(mesh,doflocs=prev["geometry"][:])
        v = prev["velocity"][:]; t = float(prev.attrs["time_s"])
        coating = Coating.from_dict(json.loads(prev.attrs["coating"]))
        prior = json.loads(prev.attrs["diagnostics"])
        dt_next = min(c.dt,prior.get("dt_s",c.dt) or c.dt)
        heartbeat(SimpleNamespace(step=n,time=t),str(native))
        started = time.monotonic()
        reference=fb.initial_mesh(c)
        while t < c.t_end-1e-13 and (stop_after is None or n<stop_after):
            # Land on physical output times even after an adaptive rejection.
            output_time=(np.floor((t+1e-11)/.01)+1)*.01
            dt = min(dt_next,c.t_end-t,output_time-t)
            start_mesh,start_v,rezoning=mesh,v,None
            last_rezone=prior.get("last_rezone_time_s",0.)
            if c.rezone_interval_s>0 and t>=last_rezone+c.rezone_interval_s-1e-11:
                from .free_boundary_remesh import rezone
                try:
                    start_mesh,start_v,rezoning=rezone(reference,mesh,v,c,orientation)
                except (ValueError,RuntimeError) as error:
                    bad=h["rejected_trials"].create_group(f"{len(h['rejected_trials']):08d}")
                    bad.attrs.update(from_step=n,time_s=t,error=str(error),operation="rezoning",accepted=False)
                    h.flush();heartbeat(SimpleNamespace(step=n,time=t),str(native))
                    raise RuntimeError(f"PW2 rezoning recovery required at t={t}: {error}") from error
                last_rezone=t
                print(json.dumps({"rezoned_at_s":t,"volume_change":rezoning["volume_change"],
                    "momentum_error":rezoning["momentum_error"],"energy_work":rezoning["energy_work"],
                    "min_jacobian":rezoning["min_jacobian"]}),flush=True)
            for trial in range(c.max_halvings+1):
                try:
                    mesh1,v1,p,row = fb.step(start_mesh,start_v,dt,c,orientation)
                    curve = fb.sample_interface(mesh1,ids,16)
                    if curve[0].min() < -1.-1e-12 or curve[0].max() > 1.+1e-12:
                        raise ValueError("free interface crosses a solid side wall")
                    break
                except (ValueError,RuntimeError) as error:
                    group = h["rejected_trials"].create_group(f"{len(h['rejected_trials']):08d}")
                    group.attrs.update(from_step=n,time_s=t,dt_s=dt,error=str(error),accepted=False)
                    h.flush(); rejected += 1
                    print(json.dumps({"rejected":True,"step":n,"time":t,"dt":dt,"error":str(error)}),flush=True)
                    if trial == c.max_halvings or dt*.5<c.minimum_dt:
                        state.update(status="NUMERICAL_RECOVERY_REQUIRED",last_error={"time_s":t,"error":str(error)})
                        heartbeat(SimpleNamespace(step=n,time=t),str(native))
                        raise RuntimeError(f"PW2 numerical recovery required at t={t}: {error}") from error
                    dt *= .5
            n += 1; t = min(c.t_end,t+dt)
            coating.accept(n,t,mesh1.p[1,ids[[0,-1]]])
            # Confirm a resolved local LEFT record only from actual accepted contact.
            if n>=2:
                before = h["states"][f"{n-2:08d}"]
                oldR = float(before["geometry"][1,ids[0]])
                lastR, newR = mesh.p[1,ids[0]],mesh1.p[1,ids[0]]
                threshold = 1e-6 # resolution criterion, not target peak height
                last_marker = max(m.height_m for m in coating.markers if m.side=="L")
                if oldR<lastR and newR<lastR and lastR>last_marker+threshold:
                    coating.record_peak(f"P{len(coating.markers)+1}","L",float(lastR),t-dt,f"{run_id}:{n-1}")
            for suffix in ("bulk_loss","left_wall_loss","geometry_kinetic_work","geometry_potential_work"):
                row["cumulative_"+suffix] = prior["cumulative_"+suffix]+row[suffix]
            row["last_rezone_time_s"]=last_rezone
            row["rezone_work"]=0. if rezoning is None else rezoning["energy_work"]
            row["cumulative_rezone_work"]=prior.get("cumulative_rezone_work",0.)+row["rezone_work"]
            row["continuous_energy_error"] = (row["energy"]-E0+row["cumulative_bulk_loss"]+row["cumulative_left_wall_loss"])/E0
            row["split_energy_error"] = row["continuous_energy_error"]-(row["cumulative_geometry_kinetic_work"]+row["cumulative_geometry_potential_work"]+row["cumulative_rezone_work"])/E0
            save_state(h,n,t,mesh1,v1,p,coating,row,rezoning)
            mesh,v,prior = mesh1,v1,row
            dt_next = min(c.dt,dt*1.25)
            if n%10==0 or t>=c.t_end-1e-13:
                print(json.dumps({"step":n,"time_s":t,"dt":dt,"R_L":float(mesh.p[1,ids[0]]),
                    "P2_error":float(mesh.p[1,ids[-1]]-op0.mesh.p[1,ids[-1]]),
                    "energy_error":row["continuous_energy_error"],"jacobian":row["min_jacobian"],
                    "iterations":row["nonlinear_iterations"],"wall_s":time.monotonic()-started}),flush=True)
                heartbeat(SimpleNamespace(step=n,time=t),str(native))
        result = {"run_id":run_id,"native":str(native),"time_s":t,"accepted_step":n,
            "completed":t>=c.t_end-1e-13,"target_5s":t>=5.-1e-13,"rejected_this_call":rejected,
            "diagnostics":prior,"coating":coating.to_dict(),"model_id":fb.MODEL_ID}
        write_json(run_dir/"summary.json",result)
        state.update(status="IN_PROGRESS",accepted_step=n,physical_time_s=t,last_checkpoint=str(native),
            next_action="Continue PW2 through numerical qualification and final package; not COMPLETE.")
        heartbeat(SimpleNamespace(step=n,time=t),str(native))
    return result


def free_boundary_main(root, cfg, args):
    mission = root/"mission/pinned_wetting"
    if args.command=="doctor":
        print(json.dumps(doctor(root),indent=2)); return 0
    if args.command=="status":
        print(json.dumps({**read_json(mission/"state.json"),"kernel_lock_busy":lock_busy(mission/"job.lock")},indent=2)); return 0
    limits = cfg["resources"]
    # Preserve exact previous status before PW2 first state mutation.
    historical = mission/"history/linear_corrected_before_PW2_state.json"
    if not historical.exists():
        if lock_busy(mission/"job.lock"): raise RuntimeError("Existing mission job active")
        write_json(historical,read_json(mission/"state.json"))
    if args.command in ("run","pilot","resume"):
        role = args.case_role or ("pilot" if args.command=="pilot" else "target")
        c = fb.Controls(nx=12*2**args.mesh_level,nz=24*2**args.mesh_level,
            dt=.0025*args.dt_scale,t_end=float(args.horizon if args.horizon is not None else (1. if args.command=="pilot" else 5.)),
            local_right_levels=args.right_refine,rezone_interval_s=args.rezone_interval,
            hydrostatic_split=args.hydrostatic_split,skew_divergence=args.skew_divergence,
            rezone_strategy=args.rezone_strategy)
        if args.run_id:
            saved = read_json(root/REL_RUNS/args.run_id/"identity.json")
            c = fb.Controls(**saved["controls"])
            role = args.run_id.split("-")[1]
        with job(root,"PW2 "+" ".join(__import__('sys').argv),limits) as (state,heartbeat):
            result = run_case(root,c,role,state,heartbeat,args.stop_after)
            print(json.dumps(result,indent=2))
        return 0
    if args.command in ("render","package","verify"):
        from .free_boundary_output import output_main
        with job(root,"PW2 "+args.command,limits) as (state,heartbeat):
            return output_main(root,args,state,heartbeat)
    raise ValueError("Unsupported PW2 command")
