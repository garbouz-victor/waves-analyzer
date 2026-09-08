"""Measured transfer summaries and plots; audit source/provenance are separate."""
import csv
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz
from sloshing.multiphase.benchmarks import divergence_transfer as b
from sloshing.multiphase.phase_rate_history import read_json,atomic_json,Journal
from sloshing.multiphase.full_rate_comparison import compare


def optional(path):return read_json(path) if Path(path).exists() else {}


def divergence_plots(rows):
    import matplotlib;matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder=b.ROOT/"analysis";folder.mkdir(parents=True,exist_ok=True)
    t=[r["time"] for r in rows]
    for name,keys in [("strong_divergence_L0.pdf",["strong_divergence_L2","orthogonal_divergence_L2"]),
            ("weak_continuity_L0.pdf",["weak_continuity_Riesz","projected_divergence_L2"])]:
        fig,ax=plt.subplots(figsize=(7,4))
        for k in keys:ax.plot(t,[r[k] for r in rows],label=k)
        ax.set_xscale("symlog",linthresh=b.POLICY["tiny_dt"]);ax.set_yscale("symlog",linthresh=1e-25)
        if name.startswith("weak"):ax.axhline(1e-12,color="red",linestyle="--",label="hard weak gate")
        ax.set(xlabel="time [s]",ylabel="L2 divergence");ax.legend();fig.tight_layout();fig.savefig(folder/name);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,k in zip(axes,("velocity_L2","E_kin","viscous_dissipation")):
        ax.plot([r[k] for r in rows],[r["strong_divergence_L2"] for r in rows]);ax.set(xlabel=k,ylabel="strong divergence L2")
    fig.tight_layout();fig.savefig(folder/"divergence_vs_hydrodynamics.pdf");plt.close(fig)


def generate():
    decision=read_json(b.ROOT/"divergence_audit/decision.json")
    moderate=optional(b.ROOT/"moderate/comparison.json");tiny=optional(b.ROOT/"target_tiny/comparison.json")
    pilot=optional(b.ROOT/"pilot/status.json");restart=optional(b.ROOT/"pilot/restart_check.json")
    boundary=optional(b.ROOT/"pilot/block_transition.json");full=optional(b.ROOT/"full_coupling_L0/status.json")
    stop=optional(b.ROOT/"STOP.json");complete=(b.ROOT/"full_coupling_L0/COMPLETE.json").exists()
    coupling={}
    if complete:
        prepared=read_json(b.historical.EQUILIBRIUM/"prepared.json")
        with np.load(b.historical.RESULTS/"linear_operator/inputs.npz") as z:phi_eq=z["phi_eq"].copy()
        coupling=compare(b.ROOT/"full_coupling_L0",b.ISO,b.schedules(),load_npz(b.historical.RESULTS/"linear_operator/M0.npz"),
                         phi_eq,prepared["mu_star"],b.POLICY)
        atomic_json(b.ROOT/"comparison/coupling.json",coupling)
        atomic_json(b.ROOT/"comparison/combined_error.json",{k:v for k,v in coupling["coupling"].items()
            if k in ("combined_D_proxy","max_combined_phi_proxy","rigorous_full_continuum_bound")})
    passed=bool(decision["passed"] and moderate.get("passed") and tiny.get("accepted") and
        restart.get("restart_check")=="passed" and boundary.get("passed") and pilot.get("level0_authorized") and
        complete and not stop and coupling.get("coupling",{}).get("passed"))
    verdict=("TAYLOR-HOOD DISCRETE INCOMPRESSIBILITY POLICY QUALIFIED BY PRESSURE-SPACE PROJECTION AND STRONG-DIVERGENCE MESH REFINEMENT; "
        "FULL CHNS PHASE-RATE BRIDGE QUALIFIED FOR THE MATCHED-DENSITY STIFF-BUMP CASE; "
        "HYDRODYNAMIC COUPLING IS BELOW THE QUALIFIED ISOLATED-CH TEMPORAL UNCERTAINTY; "
        "ISOLATED-CH TEMPORAL QUALIFICATION CONDITIONALLY TRANSFERS TO THIS FULL-CHNS BENCHMARK.") if passed else \
        "DIVERGENCE POLICY AUDIT QUALIFIED; FULL-CHNS TRANSFER NOT QUALIFIED."
    if complete and not passed:verdict+=" FULL CHNS PHASE-RATE EXECUTION PASS; ISOLATED-CH TRANSFER FAILS."
    sessions=b.sessions();audit_wall=read_json(b.ROOT/"policy/production_frozen.json")["audit_proof"]["audit_wall_s"]
    summary={"divergence_audit":decision,"moderate":moderate,"target_tiny":tiny or "not_run",
        "pilot":pilot or "not_run","restart":restart or "not_run","first_dt_boundary":boundary or "not_run",
        "full_L0":full or "not_run","full_L0_complete":complete,"coupling_transfer":"qualified" if passed else "not_qualified",
        "scientific_verdict":verdict,"overall_model":"MODEL NOT YET VALIDATED","stop":stop,
        "cost":{"audit_wall_s":audit_wall,"preliminary_wall_s":sum(r["wall_s"] for r in sessions if r["preliminary"]),
            "full_path_wall_s":sum(r["wall_s"] for r in sessions if r["full_path"]),
            "total_scientific_wall_s":audit_wall+sum(r["wall_s"] for r in sessions)},
        "remote_CI":"not observed for unpushed source; local tests only"}
    atomic_json(b.ROOT/"summary.json",summary)
    with (b.ROOT/"summary.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["metric","value"])
        for k in ("full_L0_complete","coupling_transfer","scientific_verdict","overall_model"):w.writerow([k,summary[k]])
        for k,v in summary["cost"].items():w.writerow([k,v])
    history=b.ROOT/"full_coupling_L0/scalar_history.jsonl"
    if history.exists():
        rows=Journal.read(history);divergence_plots(rows)
        from step3a8_report import plots
        plots(b.ROOT,rows,Journal.read(b.ISO/"level0/scalar_history.jsonl"),coupling)
    return summary


if __name__=="__main__":print(generate())
