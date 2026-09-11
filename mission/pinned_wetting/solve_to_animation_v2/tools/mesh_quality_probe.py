"""Read-only same-boundary mesh recovery trial; never an accepted trajectory."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

import h5py
import numpy as np
import yaml
from sloshing.pinned_wetting import free_boundary as fb
from sloshing.pinned_wetting import free_boundary_mesh_quality as quality
from sloshing.pinned_wetting import free_boundary_remesh as remesh
from sloshing.pinned_wetting.mission import job, sha256, write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native",type=Path)
    parser.add_argument("--step",type=int)
    parser.add_argument("--iterations",type=int,default=80)
    parser.add_argument("--transfer-step",type=float,default=0.,help="Optional candidate transfer and material step, NOT accepted native data")
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[4]
    output=root/"sloshing_visualization/output/pinned_wetting/free_boundary_fixed_conditions"
    limits=yaml.safe_load((root/"mission/pinned_wetting/CONFIG.yaml").read_text())["resources"]
    with job(root,"PW2 read-only exact-boundary mesh-quality recovery probe",limits) as (_,heartbeat):
        before=sha256(args.native)
        with h5py.File(args.native,"r") as h:
            controls=fb.Controls(**json.loads(h.attrs["numerical_identity"])["controls"])
            index=int(h.attrs["last_accepted_step"]) if args.step is None else args.step
            source=h["states"][f"{index:08d}"]
            physical_time=float(source.attrs["time_s"])
            reference=fb.initial_mesh(controls)
            mesh=replace(reference,doflocs=source["geometry"][:])
            velocity=source["velocity"][:]
            orientation=h["topology/orientation"][:]
            ids=h["topology/interface_nodes"][:]
        started=time.monotonic()
        candidate,statistics=quality.optimize_mesh(reference,mesh,orientation,args.iterations)
        elapsed=time.monotonic()-started
        minima,sign=fb.quadratic_minima(candidate)
        digest=sha256(Path(quality.__file__))
        stem=f"mesh_quality_probe_{args.native.parent.name}_{index}_{digest[:8]}"+("_transfer" if args.transfer_step else "")
        geometry_path=output/(stem+".npz")
        arrays={"geometry":candidate.p,"reference_geometry":reference.p,
                "original_geometry":mesh.p,"triangles":mesh.t,"interface_nodes":ids}
        transfer_report={}
        if args.transfer_step:
            try:
                candidate,new_velocity,evidence=remesh.transfer_mesh(mesh,candidate,velocity,controls,orientation)
                arrays.update({"transfer_"+key:value for key,value in evidence.items() if isinstance(value,np.ndarray)})
                transfer_report={key:value for key,value in evidence.items() if not isinstance(value,np.ndarray)}
                trial,v1,p1,row=fb.step(candidate,new_velocity,args.transfer_step,controls,orientation)
                arrays.update(trial_geometry=trial.p,trial_velocity=v1,trial_pressure=p1)
                transfer_report.update(trial_dt_s=args.transfer_step,trial_diagnostics=row,
                                       trial_P2_drift_m=float(np.max(abs(trial.p[:,ids[-1]]-mesh.p[:,ids[-1]]))))
            except (ValueError,RuntimeError) as error:
                transfer_report["error"]=str(error)
        np.savez_compressed(geometry_path,**arrays)
        report={"scope":"READ_ONLY_MESH_TRIAL_NOT_ACCEPTED_PHYSICAL_TRAJECTORY",
            "native":str(args.native.resolve()),"native_sha256":before,
            "native_unchanged":before==sha256(args.native),"source_step":index,"time_s":physical_time,
            "quality_source_sha256":digest,"optimizer_wall_s":elapsed,
            "wall_s":time.monotonic()-started,"statistics":statistics,
            "minimum_jacobian":float(minima.min()),"minimum_jacobian_cell":int(minima.argmin()),
            "orientation_preserved":bool(np.array_equal(sign,orientation)),
            "maximum_free_boundary_change_m":float(np.max(abs(candidate.p[:,ids]-mesh.p[:,ids]))),
            "candidate_geometry":str(geometry_path),"candidate_sha256":sha256(geometry_path),
            "velocity_transferred":"transfer_velocity" in arrays,
            "new_physical_step_solved":"trial_geometry" in arrays,"candidate_transfer":transfer_report,
            "remesh_source_sha256":sha256(Path(remesh.__file__)),"solver_source_sha256":sha256(Path(fb.__file__))}
        write_json(output/(stem+".json"),report)
        print(json.dumps(report,indent=2),flush=True)
        heartbeat()


if __name__=="__main__":
    main()
