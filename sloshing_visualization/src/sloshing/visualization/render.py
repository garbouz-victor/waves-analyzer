"""Explicit preview gates and streamed MP4 export, separate from FEM preparation."""

import hashlib
import json
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.backends.backend_pdf import PdfPages

from .panels import MainFigure, SurfaceFigure, TracerFigure


def signature(meta):
    payload = json.dumps(meta, sort_keys=True).encode()
    # Explorer implementation/report changes do not invalidate inspected MP4s.
    for path in (Path(__file__).parent/name for name in ("panels.py","render.py","scope.py")):
        payload += path.read_bytes()
    return hashlib.sha256(payload).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def frame_metrics(h, meta, item):
    i = item["index"]
    r = h["ranges"][i]
    return {**item, "eta_range_m": list(r[:2]), "display_grid_speed_range_m_per_s": list(r[2:4]),
            "bulk_display_grid_speed_max_m_per_s": float(r[4]),
            **{k:float(h[k][i]) for k in ("kinetic_energy","potential_energy","total_energy")},
            "surface_magnification": 1., "tracer_magnification": meta["particle_displacement_magnification"]}


def preview(cache, meta, output):
    output=Path(output)
    dest=output/"frames_preview"
    dest.mkdir(exist_ok=True)
    start=time.perf_counter()
    rows=[]
    with h5py.File(cache,"r") as h:
        fig=MainFigure(h,meta)
        for item in meta["key_frames"]:
            fig.update(item["index"]).savefig(dest/(item["name"]+".png"),dpi=100)
            rows.append(frame_metrics(h,meta,item))
            print(json.dumps(rows[-1]),flush=True)
        plt.close(fig.fig)
        for name, factory in (("vorticity",lambda:MainFigure(h,meta,omega=True)),
                              ("clean",lambda:MainFigure(h,meta,clean=True)),
                              ("tracers",lambda:TracerFigure(h,meta)),
                              ("surface",lambda:SurfaceFigure(h,meta))):
            fig=factory()
            fig.update(meta["key_frames"][1]["index"]).savefig(dest/(name+".png"),dpi=100)
            plt.close(fig.fig)
    write_json(output/"preview_metrics.json",rows)
    write_json(output/"static_preview_gate.json",{"status":"awaiting_inspection","signature":signature(meta),"runtime_s":time.perf_counter()-start})


def accept_gate(output, meta, kind):
    path=Path(output)/(kind+"_gate.json")
    value=json.loads(path.read_text())
    if value["signature"]!=signature(meta):
        raise RuntimeError("Preview predates renderer/data changes; regenerate and inspect")
    if kind == "video_preview":
        value["verified_video"] = probe_video(Path(output)/"preview_0_1s.mp4")
    value["status"]="accepted_after_inspection"
    write_json(path,value)


def require_gate(output,meta,kind):
    path=Path(output)/(kind+"_gate.json")
    if not path.exists():
        raise RuntimeError("Required preview gate missing: "+kind)
    value=json.loads(path.read_text())
    if value["status"]!="accepted_after_inspection" or value["signature"]!=signature(meta):
        raise RuntimeError("Regenerate / inspect / accept the "+kind+" before full rendering")


def probe_video(path):
    result=subprocess.run(["ffprobe","-v","error","-count_frames","-select_streams","v:0",
                           "-show_entries","stream=width,height,avg_frame_rate,nb_read_frames,duration","-of","json",str(path)],
                          check=True,capture_output=True,text=True)
    return json.loads(result.stdout)["streams"][0]


def video(cache,meta,output,kind,short=False):
    output=Path(output)
    require_gate(output,meta,"static_preview")
    if not short:
        require_gate(output,meta,"video_preview")
    names={"main":"bulk_flow_explained", "clean":"bulk_flow_clean", "vorticity":"bulk_vorticity",
           "surface":"free_surface_true_scale", "tracers":"tracer_model_comparison"}
    stem="preview_0_1s" if short else names[kind]
    path=output/(stem+".mp4")
    manifest=output/(stem+"_render.json")
    sig=signature(meta)
    if path.exists():
        if not manifest.exists():
            raise RuntimeError("Unfinished existing video; use a new output or inspect it explicitly: "+str(path))
        old=json.loads(manifest.read_text())
        if old.get("status")!="complete" or old.get("signature")!=sig:
            raise RuntimeError("Video is incomplete or stale: "+str(path))
        print("REUSE complete MP4",path,flush=True)
        return old
    start=time.perf_counter()
    info={"status":"running","signature":sig,"path":str(path),"kind":kind}
    write_json(manifest,info)
    try:
        with h5py.File(cache,"r") as h:
            indices=np.flatnonzero(h["times"][:] <= (1.+1e-12 if short else 5.+1e-12))
            if kind=="surface":fig=SurfaceFigure(h,meta)
            elif kind=="tracers":fig=TracerFigure(h,meta)
            else:fig=MainFigure(h,meta,clean=kind=="clean",omega=kind=="vorticity")
            writer=FFMpegWriter(fps=meta["fps"],codec="libx264",metadata={"title":stem,"comment":meta["interpretation"]},
                               extra_args=["-preset","veryfast","-crf","20","-pix_fmt","yuv420p","-threads","2","-movflags","+faststart"])
            with writer.saving(fig.fig,str(path),dpi=100):
                last=time.perf_counter()
                for n,i in enumerate(indices):
                    fig.update(i)
                    writer.grab_frame()
                    if n%100==0 or time.perf_counter()-last>25:
                        print(f"{stem}: frame {n+1}/{len(indices)}, t={h['times'][i]:.3f}, elapsed={time.perf_counter()-start:.1f}s",flush=True)
                        last=time.perf_counter()
            plt.close(fig.fig)
        check=probe_video(path)
        if int(check["nb_read_frames"])!=len(indices) or (check["width"],check["height"])!=(1920,1080):
            raise RuntimeError("Encoded video dimensions/frame count mismatch")
        info.update(status="complete",runtime_s=time.perf_counter()-start,bytes=path.stat().st_size,video=check)
        write_json(manifest,info)
        if short:
            write_json(output/"video_preview_gate.json",{"status":"awaiting_inspection","signature":sig,"video":check})
    except BaseException as exc:
        info.update(status="failed",failure=repr(exc))
        write_json(manifest,info)
        raise
    return info


def storyboard(cache,meta,output):
    output=Path(output)
    with h5py.File(cache,"r") as h, PdfPages(output/"key_phases.pdf") as pdf:
        fig=MainFigure(h,meta)
        for item in meta["key_frames"]:
            fig.update(item["index"])
            fig.phase.set_text(item["phase"]+" (nearest saved snapshot)")
            pdf.savefig(fig.fig,dpi=100)
        plt.close(fig.fig)
    # Compact contact sheet from the inspected key-phase images.
    from PIL import Image
    sheet=Image.new("RGB",(1920,1620),"white")
    for j,item in enumerate(meta["key_frames"]):
        with Image.open(output/"frames_preview"/(item["name"]+".png")) as im:
            thumb=im.resize((960,540),Image.Resampling.LANCZOS if hasattr(Image,"Resampling") else Image.LANCZOS)
            sheet.paste(thumb,((j%2)*960,(j//2)*540))
    sheet.save(output/"key_phases.png")


def collect_summary(output,meta):
    output=Path(output)
    manifests=[json.loads(p.read_text()) for p in sorted(output.glob("*_render.json"))]
    value={"interpretation":meta["interpretation"],"source_fine":meta["sources"]["fine"],
           "renders":manifests,"historical_step16_verdict":"unchanged: NOT QUALIFIED FOR PHYSICAL ANIMATION under global slope policy"}
    write_json(output/"render_summary.json",value)
    return value
