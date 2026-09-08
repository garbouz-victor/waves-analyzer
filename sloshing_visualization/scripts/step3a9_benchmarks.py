"""Bounded divergence audit. PART B is forbidden until the frozen rule passes."""
import argparse
from sloshing.multiphase.benchmarks import divergence_ch as b

p=argparse.ArgumentParser()
p.add_argument("--stage",required=True,choices=["preflight","freeze","rejected","prepare","M0","M1","M2","decision","audit"])
args=p.parse_args()
if args.stage=="preflight":b.preflight()
elif args.stage=="freeze":b.freeze()
elif args.stage=="rejected":b.rejected_audit()
elif args.stage=="prepare":
    for i in range(3):b.prepare_mesh(i)
elif args.stage in ("M0","M1","M2"):b.mesh_step(int(args.stage[1]))
elif args.stage=="decision":print(b.decision(),flush=True)
elif args.stage=="audit":
    b.freeze();b.rejected_audit()
    for i in range(3):b.prepare_mesh(i)
    for i in range(3):b.mesh_step(i)
    print(b.decision(),flush=True)
