#!/usr/bin/env python3
"""Measured STEP3A7 report; can truthfully report a failed/incomplete prefix."""
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.sparse import load_npz
from sloshing.multiphase.phase_rate_history import read_json,atomic_json,Journal
from sloshing.multiphase.phase_rate_analysis import compare_complete_levels
from sloshing.multiphase.benchmarks.phase_rate_series import ROOT,schedules

SUCCESS=("NONLINEAR ISOLATED CH STIFF TRANSIENT TEMPORALLY QUALIFIED ON THE CERTIFIED MESH; "
         "FULL CHNS TRANSFER AND REMAINING STEP 3A PHYSICAL BENCHMARKS STILL REQUIRED")


def fmt(v):
    if v is None: return "NOT RUN / N/A"
    if isinstance(v,float): return f"{v:.12g}"
    return str(v)


def table(headers, rows):
    return "| "+" | ".join(headers)+" |\n| "+" | ".join(["---"]*len(headers))+" |\n"+"\n".join(
        "| "+" | ".join(fmt(x) for x in row)+" |" for row in rows)+"\n"


def plots(root, schedules, result, pilot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output=root/"analysis"; output.mkdir(parents=True,exist_ok=True)
    histories=[]
    for i in range(3):
        p=root/f"level{i}/scalar_history.jsonl"
        histories.append(Journal.read(p) if p.exists() else [])
    def save(fig,name):
        fig.tight_layout(); fig.savefig(output/(name+".pdf")); plt.close(fig)
    for log in (False,True):
        fig,ax=plt.subplots()
        for i,rows in enumerate(histories):
            if not rows: continue
            selected=rows[1:] if log else rows
            ax.plot([r["time"] for r in selected],[r["CH_dissipation"] for r in selected],label=f"L{i} ({result['levels'][i]['execution_status']})")
        if log: ax.set_xscale("log")
        ax.set(xlabel="t [s]",ylabel="D_CH [W/m]",title="Nonlinear isolated CH: measured physical power")
        if any(histories): ax.legend()
        save(fig,"dissipation_log_time" if log else "dissipation")
    fig,axes=plt.subplots(3,1,figsize=(8,9),sharex=True)
    for i,rows in enumerate(histories):
        if not rows: continue
        t=[r["time"] for r in rows]; E0=rows[0]["E_total"]
        for ax,y in zip(axes,([r["cumulative_CH_dissipation"] for r in rows],
                             [E0-r["E_total"] for r in rows],[r["energy_budget_defect"] for r in rows])):
            ax.plot(t,y,label=f"L{i}")
    for ax,label in zip(axes,("physical integral D_CH [J/m]","-[F(t)-F(0)] [J/m]","continuum budget defect [J/m]")):
        ax.set_ylabel(label)
        if any(histories): ax.legend()
    axes[-1].set_xlabel("t [s]"); save(fig,"energy")
    fig,ax=plt.subplots()
    for i,rows in enumerate(histories):
        if rows: ax.plot([r["time"] for r in rows],[r["cells_across_transition_certified_min"] for r in rows],label=f"L{i}")
    ax.axhline(8.,color="black",ls="--",label="hard gate = 8"); ax.legend()
    ax.set(xlabel="t [s]",ylabel="certified transition cells",title="EVERY accepted state (not sparse field snapshots)")
    save(fig,"resolution")
    fig,ax=plt.subplots()
    for i,rows in enumerate(histories):
        if len(rows)>1: ax.plot([r["time"] for r in rows[1:]],[r["solver"]["snes_iterations"] for r in rows[1:]],".",label=f"L{i}")
    for end in schedules[0].block_ends[:-1]: ax.axvline(schedules[0].parents[end-1].end,color="gray",alpha=.15)
    ax.set(xscale="log",xlabel="t [s]",ylabel="Newton iterations",title="Production P0 controls, previous unscaled rate guess")
    if any(histories): ax.legend()
    save(fig,"solver_iterations")
    fig,ax=plt.subplots()
    common=result.get("common_times",[])
    for i,label in enumerate(("L0–L1","L1–L2")):
        if common: ax.plot([r["time"] for r in common],[r["pairs"][i]["phi_L2"] for r in common],label=label)
    if common: ax.legend()
    else: ax.text(.5,.5,"NOT RUN: three complete levels required",ha="center",transform=ax.transAxes)
    ax.set(xlabel="common exact t [s]",ylabel="phi difference, consistent P2 L2")
    save(fig,"phi_convergence")
    fig,ax=plt.subplots()
    vals=[pilot.get("cost_projection",{}).get("projected_level0_wall_s")]+[
        r.get("wall_s") if r["execution_status"]=="complete" else None for r in result["levels"]]
    labels=["pilot forecast L0","actual L0","actual L1","actual L2"]
    for i,v in enumerate(vals):
        if v is not None: ax.bar(i,v/3600.)
        else: ax.text(i,0.,"NOT COMPLETE",rotation=90,va="bottom",ha="center",fontsize=8)
    ax.set_xticks(range(4),labels); ax.set_ylabel("wall hours, NOT CPU hours"); save(fig,"cost")
    fig,axes=plt.subplots(3,1,figsize=(8,9),sharex=True)
    for i,rows in enumerate(histories):
        if len(rows)<2: continue
        t=[r["time"] for r in rows[1:]]
        for ax,key in zip(axes,("rate_L2","dt_times_rate_max_abs","rate_max_abs")):
            ax.plot(t,[r["solver"][key] for r in rows[1:]],label=f"L{i}")
    for ax,label in zip(axes,("rate L2 [1/s]","dt max|rate|","max|rate| [1/s]")):
        ax.set_ylabel(label); ax.set_xscale("log")
        if any(histories): ax.legend()
    axes[-1].set_xlabel("t [s]"); save(fig,"phase_rate")


def generate():
    root=ROOT; ss=schedules(); frozen=read_json(root/"policy/frozen.json")
    result=compare_complete_levels(root,ss,load_npz("validation_results/step3a3/linear_operator/M0.npz"),frozen["policy"]["inherited"])
    pilot=read_json(root/"pilot/status.json") if (root/"pilot/status.json").exists() else {}
    boundary=read_json(root/"block_transition/status.json") if (root/"block_transition/status.json").exists() else {}
    restart=read_json(root/"pilot/restart_check.json") if (root/"pilot/restart_check.json").exists() else {}
    stop=read_json(root/"STOP.json") if (root/"STOP.json").exists() else None
    prerequisites=bool(pilot.get("level0_authorized") and boundary.get("passed") and restart.get("restart_check")=="passed" and not stop)
    passed=result["overall_temporal_qualification"] and prerequisites
    verdict=SUCCESS if passed else "MULTI-STEP ISOLATED CH TEMPORAL QUALIFICATION INCOMPLETE OR FAILED"
    if sum(r["execution_status"]=="complete" for r in result["levels"])==2 and not passed:
        verdict="TWO-LEVEL EVIDENCE ONLY; TEMPORAL QUALIFICATION INCOMPLETE"
    test_results={}
    for path in sorted((root/"tests").glob("*.xml")):
        suites=list(ET.parse(path).getroot().iter("testsuite"))
        test_results[path.name]={k:sum(int(s.get(k,0)) for s in suites) for k in ("tests","failures","errors","skipped")}
    session_values=[read_json(p) for p in root.glob("**/sessions/*/timing.json")]
    result.update(pilot=pilot,block_transition=boundary,restart=restart,stop=stop,
        overall_temporal_qualification=passed,temporal_verdict=verdict,model_verdict="MODEL NOT YET VALIDATED",
        full_CHNS_transfer="not_qualified",tests=test_results,remote_CI="not run: no push",
        total_scientific_wall_s=sum(s["wall_s"] for s in session_values),
        schedules=[{"count":s.nsteps,"sha256":s.sha256,"parent_sha256":s.parent_sha256} for s in ss],
        frozen_execution_provenance=frozen)
    (root/"analysis").mkdir(exist_ok=True)
    atomic_json(root/"analysis/convergence.json",result)
    atomic_json(root/"summary.json",result)
    with (root/"summary.csv").open("w") as stream:
        keys=["level","execution_status","accepted_steps","integrated_D_CH","DeltaF","final_budget_defect","budget_over_DeltaF","min_certified_cells","max_mass_domain","wall_s"]
        writer=csv.DictWriter(stream,fieldnames=keys); writer.writeheader()
        for i,r in enumerate(result["levels"]): writer.writerow({k:(i if k=="level" else r.get(k)) for k in keys})
    convergence=result.get("convergence") if isinstance(result.get("convergence"),dict) else {}
    gaps=result.get("endpoint_phi_gaps",[None,None]); common=result.get("max_common_phi_gaps",[None,None])
    D=convergence.get("D_relative_gaps",[None,None]); F=convergence.get("energy_gaps",[None,None])
    comparison_rows=[["final phi L2",*gaps,convergence.get("gates",{}).get("endpoint_phi_converges")],
        ["max common-time phi L2",*common,"diagnostic"],["integrated D relative",*D,convergence.get("gates",{}).get("last_D_difference")],
        ["final F difference",*F,convergence.get("gates",{}).get("endpoint_energy_converges")],
        ["budget trend",None,None,convergence.get("gates",{}).get("budget_decreases")]]
    with (root/"analysis/convergence.csv").open("w") as stream:
        writer=csv.writer(stream); writer.writerow(["observable","L0-L1","L1-L2","passes"]); writer.writerows(comparison_rows)
    plots(root,ss,result,pilot)
    sections=[]
    def section(title,text): sections.append(f"## {len(sections)+1}. {title}\n\n{text}\n")
    section("Goal and inherited evidence","STEP3A6 proved one tiny BE step and absolute-state cancellation, not multi-step temporal convergence. This iteration tests nonlinear isolated CH only. Historical STEP3A2–6 files and FAIL/checkpoint verdicts are unchanged.")
    section("Immutable physics/source",f"Execution base: `{frozen['implementation']['execution_base_HEAD']}`; source `{frozen['implementation']['multiphase_source_sha256']}`; runner `{frozen['implementation']['runner_source_sha256']}`; policy `{frozen['policy_sha256']}`. These identify the dirty execution snapshot, not the later result commit. Inputs, P2 mesh, quad12, physical parameters and original production Newton controls remain fixed.\n\n"+table(["input","SHA"],list(frozen["policy"]["hashes"].items())))
    section("Level0 schedule and derived levels",table(["level","intervals","derived schedule SHA"],[[i,s.nsteps,s.sha256] for i,s in enumerate(ss)])+"\nParent plan hash is unchanged. Each parent endpoint is retained exactly; child solver dt=parent dt/factor. Clock-rounding discrepancy is recorded, not used to adapt dt.")
    section("Multi-step phase-rate runner","Unknowns are retained (r_phi,mu); phi_new=phi_old+dt*r_phi. First guess is semidiscrete rate; every subsequent guess is the previous accepted rate, unscaled. All intervals use atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS. No retry or absolute/BDF2 switch.")
    section("Transaction/checkpoint design","Private physical candidate and copied persistent Diagnostics pass material/finite/mass/certified-resolution/topology/energy/work checks before publication. True rate, state/old/older, time/index, maps, diagnostic initial/previous/cumulative and discrete work sums are archived with hashes. Journal prefixes and COMPLETE markers are checked; rejected states never restart. Scalar and timing journals are append-only; selected fields/checkpoints are atomically published.")
    pmetrics={"accepted intervals":pilot.get("accepted_steps"),"final pilot time":pilot.get("final_pilot_time"),
        "first block boundary index":ss[0].block_ends[0],"crossed boundary?":pilot.get("first_block_crossed"),
        "max Newton":pilot.get("max_Newton"),"max mass/domain":pilot.get("max_mass_domain"),
        "min certified cells":pilot.get("min_certified_cells"),"max weak work":pilot.get("max_weak_work_defect"),
        "wall time including fork/init":pilot.get("cost_projection",{}).get("spent_wall_s"),
        "sec/accepted interval":pilot.get("cost_projection",{}).get("mean_s"),
        "projected L0 wall s":pilot.get("cost_projection",{}).get("projected_level0_wall_s")}
    section("Pilot",table(["metric","result"],pmetrics.items())+"\nPilot is a continuous reusable L0 prefix, not temporal qualification. "+(str(stop) if stop else ""))
    section("First block transition",json.dumps(boundary,indent=2) if boundary else "NOT RUN. Required endpoint is first block end+5, i.e. step127; no alternate-guess control elected.")
    section("Restart equivalence",json.dumps(restart,indent=2) if restart else "Production restart fork NOT RUN. Tiny separate-process restart tests are reported below, not substituted for the production check.")
    section("Cost qualification","Committed machine limits: 7200s L0, 28800s whole isolated series. Forecast=1.25*max(mean,p95) of measured accepted phase-rate intervals, plus measured spent/setup. No historical CHNS seconds/step. Pilot reuse is charged once. "+json.dumps(pilot.get("cost_projection",{})))
    metric_keys=["execution_status","accepted_steps","min_dt","max_dt","wall_s","max_Newton","max_mass_domain","min_certified_cells","final_F","integrated_D_CH","DeltaF","final_budget_defect","budget_over_DeltaF"]
    level_table=table(["metric","L0","L1","L2"],[[key,*[r.get(key) for r in result["levels"]]] for key in metric_keys])
    section("Level0",level_table+"\nOnly COMPLETE levels enter formal convergence. No partial final horizon is silently substituted.")
    section("Level1",f"Status: {result['levels'][1]['execution_status']}. Fresh identical t=0; every L0 interval split exactly in two. See level table.")
    section("Level2",f"Status: {result['levels'][2]['execution_status']}. Fresh identical t=0; every L0 interval split exactly in four. No extrapolation substitutes for this level.")
    section("Mass and resolution over all levels","Both initial and target mass/domain are checked every accepted state, as are certified P2 cells>=8 and two wall crossings. Full minima/worst-cell witnesses are in summary.json. Figure: analysis/resolution.pdf.")
    section("D_CH trajectories","Independent physical power includes t=0 and every interval. Figures: dissipation.pdf and dissipation_log_time.pdf. Isolated CH has u=0 by construction; no measured full-CHNS zero powers are claimed.")
    section("Integrated D_CH convergence",table(["observable","L0–L1","L1–L2","passes?"],comparison_rows)+"\nPhysical integral uses trapezoids, never -DeltaF. Largest/first interval fractions are diagnostics in summary.json.")
    section("Energy-budget convergence","Budget=F(t)-F(0)+integrated D_CH. Shared positive-scale/change/symmetric normalizations are retained; finest nonzero abs(B)/abs(DeltaF)<=5% is separately required without hiding it behind a significance floor. Figure: energy.pdf.\n\n"+json.dumps(convergence.get("gates",{})))
    section("Endpoint phi convergence",f"Consistent P2 mass norm gaps: {gaps}. Both field and dissipation convergence are mandatory, not energy-only evidence.")
    section("Common-time convergence",f"Max phi gaps: {common}. Exact shared accepted endpoints; no interpolation. Phi/mu/F/D/integral comparisons are in convergence.json and phi_convergence.pdf.")
    section("Observed temporal orders",json.dumps(convergence.get("observed_orders",{}))+". Diagnostic only; no demand for exactly first order and no fitted thresholds.")
    work_rows=[]
    for i,r in enumerate(result["levels"]):
        for n,w in enumerate(r.get("first_work",[]),1): work_rows.append([f"L{i}:{n}",*[w[k] for k in ("dt","Delta_E","D_old","D_new","B_BE","weak_work_defect","trapezoid_gap","continuum_local_defect")]])
    section("BE work diagnostics",table(["interval","dt","Delta_E","Dold","Dnew","B_BE","weak defect","trap gap","continuum defect"],work_rows)+"\nSigned numerical B_BE is a discrete BE work remainder, NOT physical heat. It is never added to physical dissipation. Basic work is measured every step; expanded same-state residual audit is sparse and read-only, never an absolute-residual acceptance veto.")
    section("Solver statistics","Newton distributions, maximum/index, KSP failures, rejected count, all dt transitions and rate/guess norms are in summary.json. Figures: solver_iterations.pdf and phase_rate.pdf.\n\nLocal test records:\n\n"+table(["run","tests","failed","errors","skipped"],[[name,*[v[k] for k in ("tests","failures","errors","skipped")]] for name,v in test_results.items()])+"\nRemote GitHub CI NOT RUN (no push). Production 96x48 is excluded from CI.")
    section("Computational cost",f"Measured scientific wall time including initialization/fork: {result['total_scientific_wall_s']:.6g}s. Not CPU time. Stage/session records separate initialization/JIT, inclusive PETSc events, SNES, diagnostics, resolution, work, sparse residual audit and archival I/O. Figure: cost.pdf. Timing-record write itself is excluded from per-step sample but included in session wall time.")
    not_run=[f"level{i} complete trajectory" for i,r in enumerate(result["levels"]) if r["execution_status"]!="complete"]
    section("Remaining limitations","NOT RUN / not qualified: "+", ".join(not_run+["full-CHNS rate solver","CHNS algebraic/multi-step verification","targeted coupling probe","moving 90→60 contact line","theta90/120","epsilon/slip/mobility sensitivity","falling-film PDE","gravity/tank/water-air"])+". Isolated results do not establish full CHNS transfer.\n\n"+("STOP reason: "+json.dumps(stop) if stop else ""))
    section("Scientific verdict",verdict+".\n\n**MODEL NOT YET VALIDATED**.")
    Path("STEP3A7_REPORT.md").write_text("# STEP 3A.7 measured report\n\n"+"\n".join(sections))
    return result


if __name__=="__main__": generate()
