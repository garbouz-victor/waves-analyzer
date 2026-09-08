"""New read-only cost/transfer results; never overwrite a historical report."""
import csv
from pathlib import Path
import shutil
import numpy as np
from step3a10_admin import ROOT, OLD, TRAJECTORY, read, sha, numerical_guard, prefix_guard, lineage


def cost_plots(decision, actual=None):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out=ROOT/"analysis"; out.mkdir(parents=True,exist_ok=True)
    classes=decision["event_counts"]; labels=list(classes); x=np.arange(len(labels))
    fig,ax=plt.subplots(figsize=(8,4))
    ax.bar(x-.2,[classes[k]["observed"] for k in labels],.4,label="observed 1..127")
    ax.bar(x+.2,[classes[k]["future"] for k in labels],.4,label="future 128..1068")
    ax.set_xticks(x,labels,rotation=20);ax.set_yscale("log");ax.set_ylabel("exact scheduled event count");ax.legend()
    fig.tight_layout();fig.savefig(out/"cost_event_classes.pdf");plt.close(fig)
    tests=decision["backtests"];labels=[f"{r['validation_start']}–{r['validation_end']}" for r in tests];x=np.arange(len(labels))
    fig,ax=plt.subplots(figsize=(7,4))
    ax.bar(x-.2,[r["predicted_interval_s"] for r in tests],.4,label="predeclared upper forecast")
    ax.bar(x+.2,[r["actual_interval_s"] for r in tests],.4,label="measured interval wall")
    ax.set_xticks(x,labels);ax.set_ylabel("seconds (wall)");ax.legend();fig.tight_layout()
    fig.savefig(out/"cost_forecast_backtest.pdf");plt.close(fig)
    vals=[decision["old_forecast"]["projected_full_L0_wall_s"],decision["forecast"]["projected_complete_L0_s"]]
    labels=["old whole-step p95","schedule-aware"]
    if actual is not None: vals.append(actual);labels.append("actual complete L0")
    fig,ax=plt.subplots(figsize=(7,4));ax.bar(labels,vals);ax.axhline(14400,color="red",ls="--",label="unchanged L0 cap")
    ax.set_ylabel("seconds (wall)");ax.legend();fig.tight_layout();fig.savefig(out/"forecast_vs_actual.pdf");plt.close(fig)


def generate():
    from sloshing.multiphase.phase_rate_history import atomic_json, Journal
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    f=numerical_guard(); prefix_guard(after_continuation=True)
    decision=read(ROOT/"cost_audit/decision.json")
    controller_path=ROOT/"continuation/controller_status.json"
    controller=read(controller_path) if controller_path.exists() else {"execution_status":"not_run"}
    complete=(TRAJECTORY/"COMPLETE.json").exists()
    coupling=None; full=None; actual=None; energy=None
    if complete:
        from scipy.sparse import load_npz
        from sloshing.multiphase.full_rate_comparison import compare
        prepared=read(d.historical.EQUILIBRIUM/"prepared.json")
        with np.load(d.historical.RESULTS/"linear_operator/inputs.npz") as z: phi_eq=z["phi_eq"].copy()
        coupling=compare(TRAJECTORY,d.ISO,d.schedules(),load_npz(d.historical.RESULTS/"linear_operator/M0.npz"),phi_eq,prepared["mu_star"],d.POLICY)
        # compare() re-verifies COMPLETE, all histories and common/final field archives.
        rows=Journal.read(TRAJECTORY/"scalar_history.jsonl"); full=read(TRAJECTORY/"status.json")
        actual=d.used("full_path"); energy=coupling["energy"]
        energy.update(max_decomposition_roundoff=max(abs(r["work"]["decomposition_roundoff"]) for r in rows[1:]),
                      max_work_identity_roundoff=max(abs(r["work"]["weak_work_identity_roundoff"]) for r in rows[1:]))
        def extreme(key,mode=max,absolute=False):
            row=mode(rows,key=lambda r:abs(r[key]) if absolute else r[key])
            return {"value":abs(row[key]) if absolute else row[key],"step":row["step"],"time":row["time"]}
        full["mass_extremum"]=extreme("mass_error_relative_to_domain",absolute=True)
        full["resolution_extremum"]=extreme("cells_across_transition_certified_min",min)
        newton=max(rows[1:],key=lambda r:r["solver"]["snes_iterations"])
        full["maximum_Newton"]={"iterations":newton["solver"]["snes_iterations"],"step":newton["step"],"time":newton["time"]}
        audits=[{"step":r["step"],"time":r["time"],"values":{k:v for k,v in r.items() if "advec" in k or "diffus" in k}}
                for r in rows if r["step"] in d.schedules()[0].audit_indices]
        atomic_json(ROOT/"analysis/phase_advection_diffusion.json",audits)
        atomic_json(ROOT/"analysis/coupling.json",coupling)
        atomic_json(ROOT/"analysis/full_divergence.json",full["divergence_statistics"])
        atomic_json(ROOT/"analysis/full_energy.json",energy)
        # User asks sum of separate maxima, not only the maximum of pointwise sums.
        measured=coupling["coupling"]["measured"]
        proxies={"E_combined_phi":measured["phi_coupling"]+measured["phi_temporal"],
                 "E_combined_D":measured["D_coupling"]+measured["D_temporal"],
                 "max_pointwise_combined_phi":coupling["coupling"]["max_combined_phi_proxy"],
                 "label":"practical conservative uncertainty proxy; not a rigorous continuum bound"}
        atomic_json(ROOT/"analysis/combined_error.json",proxies)
        from step3a8_report import plots
        plots(ROOT,rows,Journal.read(d.ISO/"level0/scalar_history.jsonl"),coupling)
        shutil.copyfile(ROOT/"analysis/budget.pdf",ROOT/"analysis/energy_budget.pdf")
        import matplotlib;matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for name,keys in (("strong_divergence_L0.pdf",["strong_divergence_L2","orthogonal_divergence_L2"]),
                          ("weak_continuity_L0.pdf",["weak_continuity_Riesz","projected_divergence_L2"])):
            fig,ax=plt.subplots(figsize=(7,4))
            for key in keys:ax.plot([r["time"] for r in rows],[r[key] for r in rows],label=key)
            ax.set_xscale("symlog",linthresh=d.POLICY["tiny_dt"]);ax.set_yscale("symlog",linthresh=1e-25)
            if name.startswith("weak"):ax.axhline(d.POLICY["weak_continuity"],color="red",ls="--",label="hard weak gate")
            ax.set(xlabel="physical time [s]",ylabel="L2 divergence");ax.legend();fig.tight_layout()
            fig.savefig(ROOT/"analysis"/name);plt.close(fig)
    inherited=bool(read(OLD/"divergence_audit/decision.json")["passed"] and read(OLD/"moderate/comparison.json")["passed"] and
        read(OLD/"target_tiny/comparison.json")["accepted"] and read(OLD/"pilot/restart_check.json")["restart_check"]=="passed" and
        read(OLD/"pilot/block_transition.json")["passed"])
    passed=bool(decision["authorized"] and inherited and complete and controller["execution_status"]=="complete" and
                coupling and coupling["coupling"]["passed"] and d.used("full_path")<=14400 and
                d.used()+decision["forecast"]["calibration_administrative_wall_s"]<=19800)
    costs={"old_complete_L0_forecast_s":decision["old_forecast"]["projected_full_L0_wall_s"],
           "new_complete_L0_forecast_s":decision["forecast"]["projected_complete_L0_s"],
           "actual_complete_L0_s":actual,"actual_primary_L0_spent_s":d.used("full_path"),
           "actual_total_scientific_s":d.used()+decision["forecast"]["calibration_administrative_wall_s"],
           "forecast_over_actual":decision["forecast"]["projected_complete_L0_s"]/actual if actual else None}
    atomic_json(ROOT/"analysis/actual_vs_forecast.json",costs);cost_plots(decision,actual)
    if passed:
        verdict=("FULL CHNS PHASE-RATE CONTINUOUS L0 TRAJECTORY COMPLETE; HYDRODYNAMIC COUPLING IS BELOW THE QUALIFIED "
                 "ISOLATED-CH TEMPORAL UNCERTAINTY FOR THE MATCHED-DENSITY ZERO-FORCE STIFF-BUMP CASE; "
                 "ISOLATED-CH TEMPORAL QUALIFICATION CONDITIONALLY TRANSFERS TO THIS FULL-CHNS BENCHMARK.")
    elif not decision["authorized"]:
        verdict="SCHEDULE-AWARE COST AUDIT DID NOT AUTHORIZE COMPLETION; FULL-CHNS TRANSFER REMAINS INCOMPLETE."
    elif complete:verdict="FULL CHNS L0 EXECUTION QUALIFIED; ISOLATED-CH TRANSFER DOMINANCE FAILED."
    else:verdict="AUTHORIZED FULL CHNS CONTINUATION INCOMPLETE; FULL-CHNS TRANSFER NOT QUALIFIED."
    summary={"administrative_guard":"passed",**lineage(f),"cost_authorized":decision["authorized"],
             "controller":controller,"full_L0_complete":complete,"full_L0":full,"coupling_transfer":passed,
             "costs":costs,"scientific_verdict":verdict,"overall_model":"MODEL NOT YET VALIDATED",
             "remote_CI":"not observed; no push","old_reports_and_STOPs":"unchanged historical evidence"}
    atomic_json(ROOT/"resumption_summary.json",summary)
    with (ROOT/"resumption_summary.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["metric","value"])
        for k in ("administrative_guard","cost_authorized","full_L0_complete","coupling_transfer","scientific_verdict","overall_model"):
            writer.writerow([k,summary[k]])
        for k,v in costs.items():writer.writerow([k,v])
    return summary


if __name__=="__main__":print(generate())
