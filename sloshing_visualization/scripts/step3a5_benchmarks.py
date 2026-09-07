#!/usr/bin/env python3
"""Exactly one explicit STEP3A5 stage; no automatic downstream PDE runs."""
import argparse
from sloshing.multiphase.benchmarks import phase_rate_ch as stages

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=stages.STAGES)
    args = parser.parse_args()
    try:
        result = stages.STAGES[args.stage]()
        if isinstance(result, dict) and result.get("qualification_status") == "failed":
            raise RuntimeError("STOP: scientific gate failed; downstream PDE stages forbidden")
    except Exception as exc:
        if stages.ACTIVE_OUTPUT is not None:
            path = stages.ACTIVE_OUTPUT/"status.json"
            current = stages.read_json(path)
            if current.get("status") == "running":
                stages.write_json(path, {"status": "stopped", "qualification_status": "failed",
                    "reason": str(exc), "stage": args.stage})
        raise
