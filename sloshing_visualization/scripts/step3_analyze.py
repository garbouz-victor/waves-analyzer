"""Read-only STEP 3A result aggregation. No solver, tank or animation calls."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parents[1]/
    "validation_results/step3/cache/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_json(path):
    return json.loads(path.read_text())


def history(path):
    with path.open() as stream:
        rows=list(csv.DictReader(stream))
    return {key:np.array([float(row[key]) for row in rows]) for key in rows[0]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,default=Path("validation_results/step3"))
    args=parser.parse_args()
    base=args.results
    flat=read_json(base/"flat/summary.json")
    operators=read_json(base/"operators_checked/summary.json")
    laplace=[read_json(base/"laplace"/name/"laplace.json")
             for name in ("resolved","radius03","radius03_dt_half")]
    contact=read_json(base/"contact/theta60_resolved/contact.json")
    startup=read_json(base/"startup_energy/summary.json")
    serial=read_json(base/"mpi/one_rank/summary.json")["final"]
    parallel=read_json(base/"mpi/two_ranks_v2/summary.json")["final"]
    mpi={k:abs(serial[k]-parallel[k]) for k in
         ("E_total","phase_mass","speed_max_dof_sample","strong_divergence_L2")}
    observations=read_json(base/"contact/theta60_resolved/observations.json")
    violations=[o for o in observations if not o["interface_resolution"]["qualified_resolution"]]
    gates={"flat_interface_h_p":flat["status"]=="passed",
           "manufactured_operator_consistency":operators["status"]=="passed",
           "laplace_multiple_radii_energy":all(r["status"]=="passed" for r in laplace),
           "contact_60_degree_run":contact["status"]=="passed",
           "contact_90_120_degree_PDE_runs":False,
           "contact_epsilon_convergence":False,
           "contact_fine_startup_energy_qualified":False,
           "spreading_slip_mobility_study":False,
           "falling_film_flux_benchmark":False,
           "film_regularization_studies":False,
           "MPI_tiny_comparison":max(mpi.values())<1e-12}
    summary={"status":"MODEL NOT YET VALIDATED" if not all(gates.values()) else
             "MODEL VALIDATED FOR TANK BRIDGE RUN", "gates":gates,
        "tank_bridge_allowed":all(gates.values()),"film_demo_allowed":False,
        "base_commit":"c6dcf3faf28ff2cf9f303e5d5d8832d9ed2f7517",
        "flat_rows":flat["rows"],"operator_rows":operators["rows"],
        "operator_scope":operators["scope"],
        "laplace_measurements":[r["measurement"] for r in laplace],
        "laplace_energy_increments":[r["run"]["max_energy_increase_step"] for r in laplace],
        "contact_angle_measured_deg":contact["measured_theta_deg"],
        "contact_mass_error_max":contact["run"]["max_mass_error_relative"],
        "contact_first_resolution_violation_s":violations[0]["time"] if violations else None,
        "contact_min_cells_across_transition":min(o["interface_resolution"]["cells_across_transition_min"] for o in observations),
        "contact_budget_relative_max":contact["run"]["max_energy_budget_relative"],
        "contact_final_cumulative_CH_unqualified_J_per_m":contact["run"]["final"]["cumulative_CH_dissipation"],
        "startup_time_rows":startup["rows"],"MPI_absolute_differences":mpi,
        "limitations":["Time-evolved tests are matched-density demonstration fluids, not water-air.",
            "Unmatched-density operators have analytic consistency checks, not a production gravity validation.",
            "Contact interface leaves the fully refined band; 8-cell policy is unchanged.",
            "Instantaneous variational wetting plus incompatible initial angle gives stiff CH startup.",
            "Current time-integrated dissipation is not qualified physical heat.",
            "BDF2 does not guarantee monotone physical energy; the second-radius gate remains failed.",
            "No falling-film PDE benchmark, epsilon/slip/mobility convergence or tank release yet.",
            "No hydrodynamic wall film or molecular wetting claim is made."]}
    source=Path(__file__).resolve().parents[1]/"src/sloshing/multiphase"
    digest=hashlib.sha256()
    for path in sorted(source.rglob("*.py")):
        digest.update(str(path.relative_to(source)).encode());digest.update(path.read_bytes())
    summary["analysis_source_sha256"]=digest.hexdigest()
    summary["provenance_note"]="Per-run configs/checkpoints are retained; analysis source hash is not retrospectively asserted to be each development run's source hash."
    (base/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    with (base/"summary.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["gate","passed"]);writer.writerows(gates.items())
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    c=history(base/"contact/theta60_resolved/history.csv")
    axes[0,0].plot(c["time"],c["E_total"],label="physical energy")
    axes[0,0].set(xlabel="time [s]",ylabel="E [J/m]");axes[0,0].legend()
    axes[0,1].plot(c["time"],c["cumulative_CH_dissipation"],color="tab:red")
    axes[0,1].set(xlabel="time [s]",ylabel="UNQUALIFIED CH time integral [J/m]",
                  title="Do not interpret this curve as physical heat")
    axes[1,0].plot([o["time"] for o in observations],
                  [o["interface_resolution"]["cells_across_transition_min"] for o in observations])
    axes[1,0].axhline(8,color="k",ls="--")
    axes[1,0].set(xlabel="time [s]",ylabel="minimum cells across phi=-0.9..0.9")
    rows=startup["rows"]
    axes[1,1].loglog([r["dt"] for r in rows],[r["energy_budget_relative_max"] for r in rows],"o-")
    axes[1,1].set(xlabel="dt [s]",ylabel="max relative energy defect",
                 title="Coarse fixed-space temporal diagnostic ONLY")
    for ext in ("pdf","png"):
        fig.savefig(base/f"gate_diagnostics.{ext}",dpi=150)
    plt.close(fig)
    print(summary["status"])
    print(json.dumps(gates,indent=2))


if __name__=="__main__":
    main()
