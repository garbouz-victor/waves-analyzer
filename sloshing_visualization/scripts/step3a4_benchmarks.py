#!/usr/bin/env python3
"""Run exactly one bounded STEP3A4 scientific stage; no automatic downstream run."""
import argparse
from sloshing.multiphase.benchmarks import cost_efficient_ch as stages

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=stages.STAGES)
    args = parser.parse_args()
    try:
        stages.STAGES[args.stage]()
    except Exception as exc:
        if stages.ACTIVE_OUTPUT is not None:
            path = stages.ACTIVE_OUTPUT/"status.json"
            current = stages.read_json(path)
            if current.get("status") == "running":
                stages.write_json(path, {"status": "stopped", "qualification_status": "failed",
                    "reason": str(exc), "stage": args.stage})
        raise
