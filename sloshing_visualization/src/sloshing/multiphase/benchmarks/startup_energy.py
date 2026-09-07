"""Fixed-space time refinement of an incompatible microscopic wetting startup.

This inexpensive, deliberately coarse diagnostic does NOT qualify contact angle
or film thickness. It isolates quadrature/time error in the discrete CH transient.
"""
import csv
import json
from pathlib import Path
from mpi4py import MPI
from ..initialization import sessile_drop
from .common import run_case


def study(config,output):
    output=Path(output)
    rows=[]
    for dt in (.0004,.0002,.0001):
        c=config.changed(nx=16,nz=8,refinement_levels=0,dt=dt,t_end=.004)
        _,history,r=run_case(c,sessile_drop(c.epsilon,.3,90.,c.z_min),output/f"dt_{dt:.4f}")
        rows.append({"dt":dt,"energy_budget_relative_max":r["max_energy_budget_relative"],
            "initial_CH_dissipation":r["initial"]["CH_dissipation"],
            "final_energy":r["final"]["E_total"],"energy_increase_max":r["max_energy_increase_step"],
            "mass_error_max":r["max_mass_error_relative"],"runtime_s":r["runtime_s"]})
    decreasing=all(b["energy_budget_relative_max"]<a["energy_budget_relative_max"] for a,b in zip(rows,rows[1:]))
    result={"status":"passed" if decreasing else "failed","budget_decreases_with_dt":decreasing,
            "rows":rows,"scope":"fixed under-resolved spatial model: temporal diagnosis ONLY"}
    if MPI.COMM_WORLD.rank==0:
        (output/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
        with (output/"summary.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
        print(json.dumps(result,indent=2),flush=True)
    return result
