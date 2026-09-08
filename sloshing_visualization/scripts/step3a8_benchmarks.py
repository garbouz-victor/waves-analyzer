#!/usr/bin/env python3
"""Explicit bounded full-CHNS validation stages; no production work in pytest."""
import argparse
from sloshing.multiphase.benchmarks import full_rate_ch as run
from sloshing.multiphase.phase_rate_history import atomic_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",required=True,choices=("audit-references","tiny-bridge","freeze","moderate",
        "target-tiny","pilot","restart-check","block-transition","full-L0","analyze","all"))
    args=parser.parse_args()
    try:
        if args.stage=="all": run.all_stages()
        elif args.stage!="analyze": getattr(run,args.stage.replace("-","_"))()
    except Exception as exc:
        atomic_json(run.ROOT/"STOP.json",{"stage":args.stage,"error":str(exc),"exception":type(exc).__name__,
            "no_retry":True,"scientific_verdict":"MODEL NOT YET VALIDATED"})
        raise
    finally:
        if args.stage not in ("audit-references","tiny-bridge","freeze"):
            from step3a8_report import generate
            generate()


if __name__=="__main__": main()
