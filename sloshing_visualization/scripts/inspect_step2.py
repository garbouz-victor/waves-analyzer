#!/usr/bin/env python3
"""Read-only source checks and complete MP4 decoding; a few extracted audit frames."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sloshing.validation.qualification.dataset import dataset_fingerprint
from sloshing.validation.qualification.fem_evaluation import FEMGeometry,evaluate_map
from sloshing.validation.contact import surface_value
from sloshing.visualization.data import grid_map
from sloshing.visualization.render import probe_video,signature,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("output/step2"))
    parser.add_argument("--video",nargs="+",help="MP4 filenames to inspect together; default: all six full outputs")
    args=parser.parse_args();out=args.output
    meta=json.loads((out/"animation_metadata.json").read_text())
    for key,stored in meta["sources"].items():
        if dataset_fingerprint(stored["path"])!=stored:
            raise RuntimeError("Source fingerprint changed: "+key)
    errors={"u":0.,"w":0.,"omega":0.,"eta":0.}
    with h5py.File(meta["sources"]["fine"]["path"],"r") as source,h5py.File(out/"display_data.h5","r") as cache:
        geometry=FEMGeometry.from_hdf5(source)
        mp=grid_map(geometry,cache["x"][:],cache["z"][:])
        mo=grid_map(geometry,cache["x"][:],cache["z"][:],degree=1)
        for i in [r["index"] for r in meta["key_frames"]]+[1000]:
            for key in ("u","w","omega"):
                actual=evaluate_map(source["fields/"+key][i],mo if key=="omega" else mp)
                stored=cache[key][i].ravel()
                np.testing.assert_allclose(stored,actual,atol=1e-10 if key!="omega" else 2e-8,rtol=1e-7)
                errors[key]=max(errors[key],float(abs(stored-actual).max()))
            eta=surface_value(source["mesh/surface_x"][:],source["fields/eta"][i],cache["eta_x"][:])
            np.testing.assert_array_equal(eta,cache["eta"][i])
        arrow_s=meta["fixed_arrow_seconds"]
        tipz=cache["arrow_z"][:]+arrow_s*cache["arrow_w"][:]
        tipx=cache["arrow_x"][:]+arrow_s*cache["arrow_u"][:]
        assert tipz.max()<=0 and tipx.min()>=-1 and tipx.max()<=1
        arrow_bounds={"x_min":float(tipx.min()),"x_max":float(tipx.max()),"z_max":float(tipz.max())}
    names=args.video or [s+".mp4" for s in ("bulk_flow_explained","bulk_flow_clean","bulk_vorticity","free_surface_true_scale","tracer_model_comparison","no_slip_boundary_layer")]
    inspection=out/"video_inspection";inspection.mkdir(exist_ok=True)
    summary_path=out/"inspection_summary.json"
    result=json.loads(summary_path.read_text()) if summary_path.exists() else {"videos":{}}
    for name in names:
        path=out/name
        manifest=json.loads((out/(path.stem+"_render.json")).read_text())
        assert manifest["status"]=="complete" and manifest["signature"]==signature(meta)
        info=probe_video(path)
        count=201 if path.stem=="preview_0_1s" else 1001
        assert int(info["nb_read_frames"])==count and (info["width"],info["height"])==(1920,1080)
        subprocess.run(["ffmpeg","-v","error","-i",str(path),"-f","null","-"],check=True)
        data=subprocess.run(["ffprobe","-v","error","-select_streams","v:0","-show_entries","frame=best_effort_timestamp_time","-of","json",str(path)],check=True,capture_output=True,text=True)
        pts=np.array([float(f["best_effort_timestamp_time"]) for f in json.loads(data.stdout)["frames"]])
        np.testing.assert_allclose(np.diff(pts),1/meta["fps"],atol=1e-10,rtol=0)
        for index in (83,count-1):
            subprocess.run(["ffmpeg","-v","error","-i",str(path),"-vf",f"select=eq(n\\,{index})","-frames:v","1","-y",str(inspection/(path.stem+f"_{index}.png"))],check=True)
        result["videos"][name]={**info,"decode_passed":True,"constant_frame_interval_passed":True,
                                "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"bytes":path.stat().st_size}
        print("VERIFIED",name,info,flush=True)
    result.update(cache_fem_max_absolute_rounding_difference=errors,all_snapshot_arrow_tip_bounds=arrow_bounds,
                  sources=meta["sources"],signature=signature(meta))
    write_json(summary_path,result)


if __name__=="__main__":main()
