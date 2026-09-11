"""One read-only remap replay, not a trajectory or a new physical result."""
import argparse
import cProfile
from dataclasses import replace
import io
import json
from pathlib import Path
import pstats
import time

import h5py
import numpy as np
import yaml
from sloshing.pinned_wetting import free_boundary as fb
from sloshing.pinned_wetting import free_boundary_remesh as remesh
from sloshing.pinned_wetting.mission import job, sha256, write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native",type=Path)
    parser.add_argument("--step",type=int,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[4]
    output=root/"sloshing_visualization/output/pinned_wetting/free_boundary_fixed_conditions"
    limits=yaml.safe_load((root/"mission/pinned_wetting/CONFIG.yaml").read_text())["resources"]
    with job(root,"PW2 read-only saved conservative-remap cost diagnostic",limits) as (_,heartbeat):
        before=sha256(args.native)
        with h5py.File(args.native,"r") as h:
            controls=fb.Controls(**json.loads(h.attrs["numerical_identity"])["controls"])
            source=h["states"][f"{args.step:08d}"]
            transferred=h["states"][f"{args.step+1:08d}"]["rezoning_before_step"]
            reference=fb.initial_mesh(controls)
            mesh=replace(reference,doflocs=source["geometry"][:])
            velocity=source["velocity"][:]
            orientation=h["topology/orientation"][:]
            expected={key:transferred[key][:] for key in ("geometry","velocity","source_cells","source_reference_coordinates")}
        profiler=cProfile.Profile()
        started=time.monotonic()
        answer=profiler.runcall(remesh.rezone,reference,mesh,velocity,controls,orientation)
        elapsed=time.monotonic()-started
        stream=io.StringIO()
        pstats.Stats(profiler,stream=stream).strip_dirs().sort_stats("cumulative").print_stats(30)
        print(stream.getvalue(),flush=True)
        digest=sha256(Path(remesh.__file__))
        stem=f"remap_profile_{args.native.parent.name}_{args.step}_{digest[:8]}"
        profiler.dump_stats(str(output/(stem+".prof")))
        evidence=answer[2]
        report={"scope":"READ_ONLY_REMAP_REPLAY_NOT_TRAJECTORY","native":str(args.native.resolve()),
            "source_step":args.step,"native_sha256":before,"native_unchanged":before==sha256(args.native),
            "remesh_source_sha256":digest,"assembly_source_sha256":sha256(Path(fb.__file__)),
            "wall_s":elapsed,"profile":stream.getvalue(),
            "max_abs_differences":{key:float(np.max(abs(evidence[key]-value))) for key,value in expected.items()},
            "evidence":{key:value for key,value in evidence.items() if not isinstance(value,np.ndarray)}}
        write_json(output/(stem+".json"),report)
        print(json.dumps({key:value for key,value in report.items() if key!="profile"},indent=2),flush=True)
        heartbeat()


if __name__=="__main__":
    main()
