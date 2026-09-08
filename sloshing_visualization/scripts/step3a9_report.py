"""Read-only summaries of measured STEP3A9 audit artifacts; never fills NOT RUN with zero."""
from pathlib import Path
from sloshing.multiphase.phase_rate_history import read_json,atomic_json

ROOT=Path("validation_results/step3a9")


def summarize():
    path=ROOT/"divergence_audit/decision.json"
    decision=read_json(path) if path.exists() else {"passed":False,"execution_status":"incomplete"}
    sessions=[read_json(p) for p in ROOT.glob("divergence_audit/**/session_status.json")]
    summary={"divergence_audit":decision,"audit_wall_s":sum(r["wall_s"] for r in sessions),
        "audit_sessions":sessions,"production_policy":"not_run","moderate":"not_run",
        "target_tiny":"not_run","pilot":"not_run","restart":"not_run","first_dt_boundary":"not_run",
        "full_L0":"not_run","coupling_transfer":"not_qualified","overall_verdict":"MODEL NOT YET VALIDATED"}
    atomic_json(ROOT/"summary.json",summary)
    return summary


if __name__=="__main__":print(summarize())
