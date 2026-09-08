#!/usr/bin/env python3
"""Explicit bounded stages; production never runs from pytest or CI."""
import argparse
from sloshing.multiphase.benchmarks import phase_rate_series as series
from sloshing.multiphase.phase_rate_history import atomic_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",required=True,choices=("freeze","pilot","restart-check","block-transition", "level0","level1","level2","analyze","all"))
    args=parser.parse_args()
    try:
        if args.stage=="all": series.all_stages()
        elif args.stage.startswith("level"): series.level(int(args.stage[-1]))
        elif args.stage!="analyze": getattr(series,args.stage.replace("-","_"))()
    except Exception as exc:
        atomic_json(series.ROOT/"STOP.json",{"stage":args.stage,"error":str(exc),"exception":type(exc).__name__,
            "scientific_verdict":"MODEL NOT YET VALIDATED","no_retry":True})
        raise
    if args.stage in ("all","analyze"):
        from step3a7_report import generate
        generate()


if __name__=="__main__": main()
