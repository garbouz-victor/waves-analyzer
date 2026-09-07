"""Gated STEP 3A.1 experiments. No tank, film or water-air entry point."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase",choices=["compatible","laplace"],required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--config",type=Path)
    parser.add_argument("--scheme",choices=["be","bdf2"],default="bdf2")
    parser.add_argument("--refinement-factor",type=int,choices=[1,2,4],default=1)
    parser.add_argument("--dt",type=float)
    parser.add_argument("--radius",type=float,default=.30)
    args=parser.parse_args()
    historical=Path("validation_results/step3").resolve()
    if args.output.resolve()==historical or historical in args.output.resolve().parents:
        parser.error("Historical STEP 3 result tree is read-only in this iteration")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.benchmarks.common import run_case
    from sloshing.multiphase.initialization import sessile_drop
    from mpi4py import MPI
    import numpy as np
    if args.phase=="compatible":
        c=ModelConfig.from_json(args.config or "configs/step3a1/contact_energy_compatible.json")
        factor=args.refinement_factor
        stages=tuple({"dt":b["dt"]/factor,"t_end":b["t_end"]} for b in c.startup_stages)
        c=c.changed(time_scheme=args.scheme,dt=c.dt/factor,startup_stages=stages,
                    refinement_circle_z=c.z_min-c.refinement_circle_radius*np.cos(np.deg2rad(c.theta_equilibrium_deg)))
        if args.dt is not None:
            parser.error("Compatible case uses its declared complete schedule, not an isolated --dt override")
        solver,history,result=run_case(c,sessile_drop(c.epsilon,c.refinement_circle_radius,
                                      c.theta_equilibrium_deg,c.z_min),args.output)
        energy=result["energy_validation"]
        checks={"mass":energy["mass_ok"],"energy_growth":energy["energy_growth_ok"],
                "energy_budget":energy["energy_budget_ok"],
                "resolution":result["interface_resolved_all_times"]}
        qualification={"qualification_status":"passed" if all(checks.values()) else "failed",
            "evidence_kind":"compatible_single_history",
            "scope":"one compatible-angle energy history; temporal/BE-vs-BDF2 series still required",
            "checks":checks,"energy_validation":energy,"run":result,
            "initialization":"direct tanh geometry with theta_initial=theta_e=60; NOT a prepared equilibrium"}
        if MPI.COMM_WORLD.rank==0:
            (args.output/"qualification.json").write_text(json.dumps(qualification,indent=2,allow_nan=False)+"\n")
            print(json.dumps({"checks":checks,"energy_validation":energy},indent=2),flush=True)
        if qualification["qualification_status"]!="passed":
            raise RuntimeError("Compatible-angle baseline FAILED: stop before Laplace/incompatible startup")
    else:
        from sloshing.multiphase.qualification_summary import read_qualification
        prerequisite=read_qualification("validation_results/step3a1/energy_compatible/qualification.json")
        if prerequisite["qualification_status"]!="passed" or prerequisite.get("evidence_kind")!="compatible_energy_series":
            raise RuntimeError("Compatible energy series gate is not passed; Laplace run blocked")
        from sloshing.multiphase.benchmarks.laplace_pressure import run
        c=ModelConfig.from_json(args.config or "configs/step3/benchmark_laplace.json")
        if args.dt is None:
            parser.error("Laplace refinement requires an explicit --dt")
        _,result=run(c.changed(dt=args.dt),args.output,args.radius)
        if result["qualification_status"]!="passed":
            raise RuntimeError("Laplace repaired gate failed; inspect measurements")


if __name__=="__main__":
    main()
