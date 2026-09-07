"""Read-only scientific plots/tables of accepted histories, including failed gates."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import simpson
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"src"))
os.environ.setdefault("MPLCONFIGDIR",str(root/"validation_results/step3a1/cache/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sloshing.multiphase.energy_validation import read_history,validate_energy
from sloshing.multiphase.settling import validate_settling


def load_json(path):
    return json.loads(path.read_text())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,default=root/"validation_results/step3a1")
    args=parser.parse_args();base=args.results
    old_folder=root/"validation_results/step3/contact/theta60_resolved"
    old=read_history(old_folder/"history.csv")
    observations=load_json(old_folder/"observations.json")
    cases=[("historical 90 to 60 (old source)",old,None)]
    for path in sorted((base/"energy_compatible").glob("*/history.csv")):
        summary=load_json(path.parent/"summary.json")
        cases.append((f"compatible 60 to 60: {path.parent.name}",read_history(path),summary))
    if len(cases)==1:
        raise ValueError("No new accepted histories to plot yet")
    table=[];quadratures=[]
    for name,rows,run in cases:
        e=validate_energy(rows);ch=e["dissipation_components"]["CH"]
        blocks=[{"offset":0,"count":len(rows)-1}] if run is None else run["provenance"]["dt_schedule"]["blocks"]
        comparison={"case":name,"scope":"offline uniform-block Simpson diagnostic; primary gate remains trapezoid"}
        for component in ("CH","viscous","slip"):
            higher=0.
            for block in blocks:
                subset=rows[block["offset"]:block["offset"]+block["count"]+1]
                times=np.array([r["time"] for r in subset])
                if len(times)<3 or not np.allclose(np.diff(times),times[1]-times[0],rtol=1e-10,atol=1e-15):
                    raise ValueError("Block Simpson requires a uniformly sampled block")
                higher+=float(simpson(np.array([r[component+"_dissipation"] for r in subset]),x=times))
            primary=e["dissipation_components"][component]["integral_J_per_m"]
            comparison[component]={"trapezoid_J_per_m":primary,"block_simpson_J_per_m":higher,
                                  "difference_J_per_m":primary-higher}
        quadratures.append(comparison)
        table.append({"case":name,"startup_dt_s":rows[1]["time"]-rows[0]["time"],
            "mesh_unknowns":82769 if run is None else run["mesh"]["total_unknowns"],
            "min_cells":min(o["interface_resolution"]["cells_across_transition_min"] for o in observations)
                if run is None else min(r["interface_cells_min"] for r in rows),
            "resolution_definition":"old centroid" if run is None else "Bernstein + seven gradient points",
            "D_CH_initial_W_per_m":ch["initial_power_W_per_m"],
            "first_interval_CH_J_per_m":ch["first_interval_J_per_m"],
            "total_CH_J_per_m":ch["integral_J_per_m"],
            "delta_E_J_per_m":e["energy_change_final_J_per_m"],
            "budget_defect_max_J_per_m":e["absolute_budget_defect"],
            "budget_defect_final_J_per_m":e["final_budget_defect_J_per_m"],
            "relative_defect_initial":e["relative_budget_defect_to_initial_scale"],
            "relative_defect_energy_change":e["relative_budget_defect_to_energy_change"],
            "CH_largest_interval_fraction":ch["largest_interval_fraction"]})
    with (base/"startup_table.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=table[0]);writer.writeheader();writer.writerows(table)
    (base/"quadrature_comparison.json").write_text(json.dumps(quadratures,indent=2,allow_nan=False)+"\n")
    for logtime in (False,True):
        fig,axes=plt.subplots(2,len(cases),figsize=(6*len(cases),7),squeeze=False,constrained_layout=True)
        for i,(name,rows,_) in enumerate(cases):
            times=np.array([r["time"] for r in rows])
            mask=times>0 if logtime else np.ones(times.size,dtype=bool)
            for component in ("CH","viscous","slip"):
                power=np.array([r[component+"_dissipation"] for r in rows])
                axes[0,i].plot(times[mask],power[mask],label=component)
            defect=np.array([r["energy_budget_defect"] for r in rows])
            axes[1,i].plot(times[mask],defect[mask],color="tab:red",label="cumulative budget defect")
            axes[0,i].set(title=name,ylabel="power [W/m]",xlabel="time [s]")
            axes[1,i].set(ylabel="budget defect [J/m]",xlabel="time [s]")
            axes[1,i].ticklabel_format(axis="y",useOffset=False)
            energy_gate=validate_energy(rows)
            change_ratio=energy_gate["relative_budget_defect_to_energy_change"]
            change_label="N/A" if change_ratio is None else f"{100*change_ratio:.4g}%"
            axes[1,i].set_title(
                f"Energy gate: {'PASS' if energy_gate['qualified'] else 'FAIL'}; "
                f"max defect / actual change: {change_label}",fontsize=9)
            axes[0,i].legend();axes[1,i].legend()
            axes[0,i].text(.03,.97,f"D_CH(0) = {rows[0]['CH_dissipation']:.6g} W/m",
                           transform=axes[0,i].transAxes,va="top",fontsize=9)
            if logtime:
                axes[0,i].set_xscale("log");axes[1,i].set_xscale("log")
        if logtime:
            fig.suptitle("t=0 cannot lie on a log axis; its power is printed and retained in every integral")
        filename="startup_power_logtime" if logtime else "startup_power_linear"
        for ext in ("pdf","png"):
            fig.savefig(base/f"{filename}.{ext}",dpi=140)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    axes[0].plot([o["time"] for o in observations],
                 [o["interface_resolution"]["cells_across_transition_min"] for o in observations],
                 label="historical moving contact: centroid indicator")
    final_audit=base/"mesh_resolution_tests/historical_final.json"
    if final_audit.exists():
        audit=load_json(final_audit)
        # This is a read-only checkpoint measurement, not a reconstructed history.
        measured=audit["new_Bernstein_indicator"]
        axes[0].plot(measured["time"],measured["cells_across_transition_min"],"rx",
                     label="old final checkpoint: new indicator (single point)")
    for name,rows,run in cases[1:]:
        axes[0].plot([r["time"] for r in rows],[r["interface_cells_min"] for r in rows],label=name)
        worst=run["worst_interface_resolution"]
        points=np.array(worst["worst_cell"]["vertices"])
        closed=np.vstack((points,points[0]))
        axes[1].plot(closed[:,0],closed[:,1],label=f"{name}, t={worst['time']:.3g} s")
    axes[0].axhline(8,color="k",ls="--",label="unchanged threshold = 8")
    axes[0].set(xlabel="time [s]",ylabel="minimum nominal transition cells")
    axes[1].set(xlabel="x [m]",ylabel="z [m]",title="New history: worst active triangle",aspect="equal")
    for ax in axes: ax.legend(fontsize=7)
    fig.suptitle("Different cases/indicator versions — not a moving-contact convergence series")
    for ext in ("pdf","png"):
        fig.savefig(base/f"resolution_history.{ext}",dpi=140)
    plt.close(fig)
    settling=validate_settling(observations,old,validate_energy(old)["E_scale_initial_J_per_m"])
    (base/"historical_settling_audit.json").write_text(json.dumps(settling,indent=2,allow_nan=False)+"\n")
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    t=np.array([o["time"] for o in observations])
    theta=np.array([o["fits"][1]["theta_deg"] for o in observations])
    crosses=np.array([sorted(p["coordinate"] for p in o["crossings"]) for o in observations])
    axes[0].plot(t,theta);axes[0].set(xlabel="time [s]",ylabel="theta [degree]",
        title=f"late dtheta/dt = {settling['angle_rate_deg_per_s']:.6g} degree/s")
    axes[1].plot(t,crosses[:,0],label="left");axes[1].plot(t,crosses[:,1],label="right")
    axes[1].set(xlabel="time [s]",ylabel="contact crossing x [m]",
        title=f"late speeds: {settling['contact_left_speed_m_per_s']:.4g}, {settling['contact_right_speed_m_per_s']:.4g} m/s")
    for ax in axes:
        ax.axvspan(settling["window_start_s"],settling["window_end_s"],color="gray",alpha=.15)
    axes[1].legend();fig.suptitle("Historical run re-audit only; no new settled moving-contact solution is claimed")
    for ext in ("pdf","png"):
        fig.savefig(base/f"settling_history.{ext}",dpi=140)
    plt.close(fig)
    print(json.dumps(table,indent=2))


if __name__=="__main__":
    main()
