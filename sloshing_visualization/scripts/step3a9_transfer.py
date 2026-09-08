"""PART B only, after a qualified audit; no automatic retries or partial transfer."""
import argparse
from sloshing.multiphase.benchmarks import divergence_transfer as b
from sloshing.multiphase.phase_rate_history import atomic_json

p=argparse.ArgumentParser()
p.add_argument("--stage",required=True,choices=["freeze","moderate","target-tiny","pilot","restart",
    "restart-fork","boundary","full-L0","analyze","all"])
args=p.parse_args()
try:
    if (b.ROOT/"STOP.json").exists():raise RuntimeError("Recorded scientific STOP forbids resumption")
    if args.stage=="freeze":b.freeze()
    elif args.stage=="moderate":b.moderate()
    elif args.stage=="target-tiny":b.target_tiny()
    elif args.stage=="pilot":b.pilot()
    elif args.stage=="restart":b.restart_check()
    elif args.stage=="restart-fork":b.execute(28,folder=b.ROOT/"pilot/restart_fork",resume=True,full_path=False)
    elif args.stage=="boundary":b.block_transition()
    elif args.stage=="full-L0":b.full_L0()
    elif args.stage=="all":b.all_stages()
    if args.stage in ("all","analyze"):
        from step3a9_transfer_report import generate
        generate()
except Exception as exc:
    if not (b.ROOT/"STOP.json").exists():
        atomic_json(b.ROOT/"STOP.json",{"stage":args.stage,"exception":type(exc).__name__,"error":str(exc),
            "implementation":b.implementation(),"overall_model":"MODEL NOT YET VALIDATED"})
    raise
