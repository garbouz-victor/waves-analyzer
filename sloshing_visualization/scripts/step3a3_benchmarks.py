"""Staged STEP 3A.3 entry point. Each stage writes only the new result tree."""
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("operator", "spectrum", "restarted-reference", "rational-reference", "low-mode", "planner"))
    args = parser.parse_args()
    from sloshing.multiphase.benchmarks import spectral_ch as stages
    actions = {"operator": stages.operator_stage, "spectrum": stages.spectrum_stage,
        "restarted-reference": stages.restarted_reference_stage,
        "rational-reference": stages.rational_reference_stage, "low-mode": stages.low_mode_stage,
        "planner": stages.planner_stage}
    try:
        result = actions[args.stage]()
    except BaseException as error:
        if stages.ACTIVE_OUTPUT is not None:
            stages.write_json(stages.ACTIVE_OUTPUT/"status.json", {
                "status": "interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                "qualification_status": "failed", "stage": args.stage,
                "error": str(error), "error_type": type(error).__name__,
                "scientific_verdict": "MODEL NOT YET VALIDATED"})
        raise
    if result["qualification_status"] != "passed":
        raise SystemExit("STOP: scientific prerequisite failed")


if __name__ == "__main__":
    main()
