"""Flat-interface h/p study from an analytic tanh, NOT discrete equilibrium."""
import csv
import json
from pathlib import Path
from mpi4py import MPI
from .common import run_case
from ..initialization import flat_interface
from ..validation_policy import MIN_TRANSITION_CELLS, EnergyPolicy


def study(config, output):
    output=Path(output)
    rows=[]
    for degree,n in ((2,16),(2,32),(2,48),(1,48)):
        c=config.changed(nx=n,nz=n,phase_degree=degree)
        _,_,r=run_case(c,flat_interface(c.epsilon),output/f"p{degree}_n{n}")
        rows.append({"phase_degree":degree,"n":n,"h":1/n,
            "spurious_speed_max":r["max_speed_m_per_s"],
            "surface_energy_error_abs":abs(r["final"]["E_interface"]-c.sigma*(c.x_max-c.x_min)),
            "mass_error_relative":r["max_mass_error_relative"],
            "energy_increase_max":r["max_energy_increase_step"],
            "energy_budget_relative":r["max_energy_budget_relative"],
            "energy_budget_qualified":r["energy_validation"]["qualified"],
            "transition_cells":r["final_interface_resolution"]["cells_across_transition_min"],
            "runtime_s":r["runtime_s"]})
    # Declared before looking at study results: monotone spurious-current and
    # interfacial-energy error reduction, mass tolerance, resolved finest mesh.
    p2=rows[:3]
    checks={"spurious_currents_decrease":all(b["spurious_speed_max"] < a["spurious_speed_max"]
        for a,b in zip(p2,p2[1:])),
        "surface_energy_error_decreases":all(b["surface_energy_error_abs"] < a["surface_energy_error_abs"]
        for a,b in zip(p2,p2[1:])),
        "mass_conservation":max(r["mass_error_relative"] for r in rows)<=EnergyPolicy().mass_relative_tolerance,
        "no_unexplained_energy_growth":max(r["energy_increase_max"] for r in rows)<=EnergyPolicy().growth_absolute_tolerance,
        "energy_budget_closure":all(r["energy_budget_qualified"] for r in rows),
        "finest_interface_resolved":p2[-1]["transition_cells"]>=MIN_TRANSITION_CELLS}
    status="passed" if all(checks.values()) else "failed"
    result={"status":status,"qualification_status":status,"checks":checks,"rows":rows,
            "scope":"flat matched-density benchmark only; not whole-model qualification"}
    if MPI.COMM_WORLD.rank==0:
        with (output/"summary.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
        (output/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(9,3.6),constrained_layout=True)
        for degree in (1,2):
            data=[r for r in rows if r["phase_degree"]==degree]
            axes[0].loglog([r["h"] for r in data],[r["spurious_speed_max"] for r in data],"o-",label=f"CG{degree} phase")
            axes[1].loglog([r["h"] for r in data],[r["surface_energy_error_abs"] for r in data],"o-")
        axes[0].set(xlabel="h [m]",ylabel="max transient spurious speed [m/s]")
        axes[1].set(xlabel="h [m]",ylabel="surface-energy error [J/m]")
        axes[0].legend()
        for ext in ("pdf","png"):
            fig.savefig(output/f"flat_convergence.{ext}",dpi=160)
        plt.close(fig)
        print(json.dumps(result,indent=2),flush=True)
    if not all(checks.values()):
        raise RuntimeError("Flat-interface gate failed; investigate before Laplace/contact benchmarks")
    return result
