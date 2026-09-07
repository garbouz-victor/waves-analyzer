"""STEP 3A benchmark CLI. Never starts a tank or visualization pipeline."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sloshing.multiphase.config import ModelConfig


def parse_arguments(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark",choices=["flat","laplace","contact","operators","startup-energy"],default="flat")
    parser.add_argument("--config",type=Path)
    parser.add_argument("--output",type=Path,default=Path("validation_results/step3/flat"))
    parser.add_argument("--nx",type=int)
    parser.add_argument("--nz",type=int)
    parser.add_argument("--dt",type=float)
    parser.add_argument("--t-end",type=float)
    parser.add_argument("--phase-degree",type=int)
    parser.add_argument("--radius",type=float)
    parser.add_argument("--theta",type=float)
    parser.add_argument("--initial-theta",type=float,default=90.)
    parser.add_argument("--epsilon",type=float)
    parser.add_argument("--mobility",type=float)
    parser.add_argument("--slip-length",type=float,dest="slip_length_m")
    parser.add_argument("--refinement-levels",type=int)
    parser.add_argument("--refinement-band",type=float,dest="refinement_band_m")
    parser.add_argument("--max-unknowns",type=int)
    execution=parser.add_mutually_exclusive_group()
    execution.add_argument("--estimate-only",action="store_true")
    parser.add_argument("--study",action="store_true",help="Run declared spatial/phase-order refinement")
    execution.add_argument("--analyze-only",action="store_true",help="Analyze complete saved Laplace fields; no PDE timesteps")
    args=parser.parse_args(argv)
    if args.analyze_only and args.benchmark!="laplace":
        parser.error("--analyze-only is implemented only for laplace; no solver was started")
    if args.study and args.benchmark not in ("flat","operators","startup-energy"):
        parser.error("--study is not implemented for this benchmark; no solver was started")
    if args.analyze_only and any(getattr(args,key) is not None for key in
            ("config","nx","nz","dt","t_end","phase_degree","radius","theta",
             "epsilon","mobility","slip_length_m","refinement_levels",
             "refinement_band_m","max_unknowns")):
        parser.error("--analyze-only uses the saved config; parameter overrides are not allowed")
    return args


def main():
    args=parse_arguments()
    default_name="flat" if args.benchmark=="operators" else args.benchmark
    if args.benchmark=="startup-energy":
        default_name="contact"
    config=ModelConfig.from_json(args.config or Path(f"configs/step3/benchmark_{default_name}.json"))
    changes={key:getattr(args,key) for key in ("nx","nz","dt","t_end","phase_degree",
             "epsilon","mobility","slip_length_m","refinement_levels","refinement_band_m","max_unknowns")
             if getattr(args,key) is not None}
    config=config.changed(**changes)
    if args.theta is not None:
        config=config.changed(theta_equilibrium_deg=args.theta)
    from sloshing.multiphase.mesh import estimate_unknowns
    if args.estimate_only:
        print(estimate_unknowns(config));return
    from sloshing.multiphase.initialization import flat_interface
    from sloshing.multiphase.benchmarks.common import run_case
    if args.benchmark=="startup-energy":
        from sloshing.multiphase.benchmarks.startup_energy import study
        result=study(config,args.output)
        if result["status"]!="passed":
            raise RuntimeError("Startup time-refinement diagnostic failed; inspect summary.json")
        return
    if args.benchmark=="operators":
        from sloshing.multiphase.benchmarks.manufactured_operators import study
        study(config,args.output)
        return
    if args.benchmark=="contact":
        from sloshing.multiphase.benchmarks.contact_angle import run
        _,result=run(config,args.output,args.initial_theta)
        if result["status"]!="passed":
            raise RuntimeError("Contact-angle gate not passed; inspect contact.json before continuing")
        return
    if args.benchmark=="laplace":
        from sloshing.multiphase.benchmarks.laplace_pressure import run,analyze_complete
        if args.analyze_only:
            result=analyze_complete(args.output)
        else:
            _,result=run(config,args.output,args.radius)
        if result["status"]!="passed":
            raise RuntimeError("Laplace gate not passed; inspect laplace.json before continuing")
        return
    if args.study:
        from sloshing.multiphase.benchmarks.static_interface import study
        study(config,args.output)
        return
    run_case(config,flat_interface(config.epsilon),args.output)


if __name__=="__main__":
    main()
