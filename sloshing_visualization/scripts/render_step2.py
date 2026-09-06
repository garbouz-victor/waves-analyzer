#!/usr/bin/env python3
"""Render existing qualified bulk fields. Missing data NEVER launches simulation."""

import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))

from sloshing.visualization.data import prepare_data
from sloshing.visualization.render import preview, video, accept_gate, collect_summary
from sloshing.visualization.storyboard import storyboard


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset",type=Path,default=Path("validation_results/animation_qualification/runs/fine-full.h5"))
    p.add_argument("--medium",type=Path,default=Path("validation_results/animation_qualification/runs/medium-full.h5"))
    p.add_argument("--output",type=Path,default=Path("output/step2"))
    for name in ("prepare","preview","preview-video","accept-preview","accept-preview-video","main","clean","vorticity","surface","tracers","storyboard","explorer","no-slip","profiles","all"):
        p.add_argument("--"+name,action="store_true")
    args=p.parse_args()
    cache,meta=prepare_data(args.dataset,args.medium,args.output)
    if args.preview:preview(cache,meta,args.output)
    if args.accept_preview:accept_gate(args.output,meta,"static_preview")
    if args.preview_video:video(cache,meta,args.output,"main",short=True)
    if args.accept_preview_video:accept_gate(args.output,meta,"video_preview")
    for kind in ("main","clean","vorticity","surface","tracers","no_slip"):
        if getattr(args,kind) or args.all:
            video(cache,meta,args.output,kind)
    if args.storyboard or args.all:storyboard(cache,meta,args.output)
    if args.profiles or args.all:
        from sloshing.visualization.no_slip_figures import static_profiles
        static_profiles(cache,meta,args.output)
    if args.explorer or args.all:
        from sloshing.visualization.explorer import explorer
        explorer(cache,meta,args.output)
    collect_summary(args.output,meta)


if __name__=="__main__":main()
