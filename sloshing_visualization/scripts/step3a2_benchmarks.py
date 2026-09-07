"""STEP 3A.2 stages; deliberately no moving-angle quench or tank command."""
import argparse
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.benchmarks.prepared_equilibrium import prepare,transient


def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="command",required=True)
    p=sub.add_parser("prepare")
    p.add_argument("--config",default="configs/step3a2/equilibrium60_symmetric.json")
    p.add_argument("--output",default="validation_results/step3a2/equilibrium60/prepared_symmetric")
    p.add_argument("--nx",type=int)
    p.add_argument("--nz",type=int)
    for command in ("preserve","perturb"):
        p=sub.add_parser(command)
        p.add_argument("--equilibrium",default="validation_results/step3a2/equilibrium60/prepared_symmetric")
        p.add_argument("--output",required=True)
        p.add_argument("--scheme",choices=("be","bdf2"),default="be")
        p.add_argument("--dt",type=float,default=1e-5)
        p.add_argument("--t-end",type=float,default=1e-3 if command=="preserve" else 1e-4)
    args=parser.parse_args()
    if args.command=="prepare":
        c=ModelConfig.from_json(args.config)
        if args.nx is not None or args.nz is not None:
            if args.nx is None or args.nz is None:
                parser.error("Specify both nx and nz")
            c=c.changed(nx=args.nx,nz=args.nz)
        result=prepare(c,args.output)
    else:
        result=transient(args.equilibrium,args.output,args.scheme,args.dt,args.t_end,args.command=="perturb")
    if result["qualification_status"]!="passed":
        raise SystemExit(2)


if __name__=="__main__":
    main()
