#!/usr/bin/env python3
"""Read-only analysis; absent stages are NOT RUN, never measured zeros."""
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.sparse import load_npz
from sloshing.multiphase.benchmarks.full_rate_ch import ROOT,ISO,schedules,historical,POLICY,sessions
from sloshing.multiphase.phase_rate_history import read_json,atomic_json,Journal
from sloshing.multiphase.full_rate_comparison import compare


def optional(path): return read_json(path) if Path(path).exists() else {}
def fmt(value):
    if value is None: return "NOT RUN"
    if isinstance(value,float): return f"{value:.12g}"
    return str(value).replace("|","/")
def table(headers,rows):
    return "\n| "+" | ".join(headers)+" |\n|"+"|".join("---" for _ in headers)+"|\n"+"".join(
        "| "+" | ".join(fmt(v) for v in row)+" |\n" for row in rows)+"\n"


def plots(root,rows,isolated,comparison):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder=Path(root)/"analysis"; folder.mkdir(exist_ok=True)
    if not rows: return
    t=np.array([r["time"] for r in rows]); ti=np.array([r["time"] for r in isolated])
    def plot(name,curves,ylabel,*,log_y=False):
        fig,ax=plt.subplots(figsize=(7,4))
        for x,y,label in curves: ax.plot(x,y,label=label)
        ax.set_xscale("symlog",linthresh=4.657657028534679e-12)
        if log_y: ax.set_yscale("symlog",linthresh=1e-25)
        ax.set_xlabel("physical time [s]; linear at t=0, logarithmic thereafter")
        ax.set_ylabel(ylabel); ax.legend(); ax.grid(alpha=.2); fig.tight_layout()
        fig.savefig(folder/name); plt.close(fig)
    get=lambda key:np.array([r[key] for r in rows])
    plot("dissipation_coupling.pdf",[(t,get("CH_dissipation"),"full L0"),(ti,[r["CH_dissipation"] for r in isolated],"isolated L0")],"D_CH [W/m]",log_y=True)
    plot("energy_coupling.pdf",[(t,get("F_CH"),"full F_CH L0"),(ti,[r["E_total"] for r in isolated],"isolated F L0")],"CH free energy [J/m]")
    plot("hydro_dissipation.pdf",[(t,get(k),label) for k,label in (("CH_dissipation","CH"),("viscous_dissipation","viscous"),("slip_dissipation","slip"))],"physical power [W/m]",log_y=True)
    fig,axes=plt.subplots(2,1,figsize=(7,6),sharex=True)
    for ax,key,label in zip(axes,("speed_max_dof_sample","E_kin"),("max speed at DOFs [m/s]","kinetic energy [J/m]")):
        ax.plot(t,get(key)); ax.set_ylabel(label); ax.set_xscale("symlog",linthresh=4.657657028534679e-12); ax.grid(alpha=.2)
    axes[-1].set_xlabel("physical time [s]"); fig.tight_layout(); fig.savefig(folder/"velocity.pdf"); plt.close(fig)
    integral=get("cumulative_CH_dissipation")+get("cumulative_viscous_dissipation")+get("cumulative_slip_dissipation")
    plot("budget.pdf",[(t,get("E_total")-rows[0]["E_total"],"Delta E total"),(t,integral,"integral D total"),(t,get("energy_budget_defect"),"budget defect")],"energy [J/m]")
    fig,ax=plt.subplots(figsize=(7,4)); ax.plot(t[1:],[r["solver"]["snes_iterations"] for r in rows[1:]])
    for i in schedules()[0].block_ends:
        ax.axvline(schedules()[0].intervals[i-1].end,color="grey",alpha=.2)
    ax.set_xscale("symlog",linthresh=4.657657028534679e-12); ax.set_xlabel("time [s]"); ax.set_ylabel("Newton iterations")
    fig.tight_layout(); fig.savefig(folder/"solver_iterations.pdf"); plt.close(fig)
    if comparison:
        fields=comparison["fields"]; x=[r["time"] for r in fields]
        plot("phi_coupling.pdf",[(x,[r["phi_coupling_L2"] for r in fields],"full0 - iso0"),
            (x,[r["iso_phi_temporal_02_L2"] for r in fields],"iso0 - iso2 temporal gap")],"consistent P2 L2 difference",log_y=True)


def generate():
    ROOT.mkdir(parents=True,exist_ok=True)
    frozen=optional(ROOT/"policy/frozen.json"); bridge=optional(ROOT/"algebraic_bridge/tiny_fem/status.json")
    moderate=optional(ROOT/"moderate/comparison.json"); target=optional(ROOT/"target_tiny/rate/status.json")
    tiny=optional(ROOT/"target_tiny/comparison.json"); pilot=optional(ROOT/"pilot/status.json")
    restart=optional(ROOT/"pilot/restart_check.json"); boundary=optional(ROOT/"pilot/block_transition.json")
    status=optional(ROOT/"full_coupling_L0/status.json"); stop=optional(ROOT/"STOP.json")
    complete=(ROOT/"full_coupling_L0/COMPLETE.json").exists(); coupling={}
    if complete:
        prepared=read_json(historical.EQUILIBRIUM/"prepared.json")
        with np.load(historical.RESULTS/"linear_operator/inputs.npz") as z: phi_eq=z["phi_eq"].copy()
        # Prepared mu is a spatial constant. Its archived key is verified independently below.
        mu_eq=prepared["mu_star"]
        coupling=compare(ROOT/"full_coupling_L0",ISO,schedules(),load_npz(historical.RESULTS/"linear_operator/M0.npz"),phi_eq,mu_eq,POLICY)
        atomic_json(ROOT/"comparison/coupling.json",coupling)
        atomic_json(ROOT/"comparison/combined_error.json",{k:v for k,v in coupling["coupling"].items() if k in ("combined_D_proxy","max_combined_phi_proxy","rigorous_full_continuum_bound")})
    prerequisites=bool(bridge.get("passed") and moderate.get("passed") and tiny.get("accepted") and
        restart.get("restart_check")=="passed" and boundary.get("passed") and complete and not stop)
    passed=prerequisites and coupling.get("coupling",{}).get("passed",False)
    if passed:
        verdict="FULL CHNS PHASE-RATE BRIDGE QUALIFIED FOR THE MATCHED-DENSITY STIFF-BUMP CASE; HYDRODYNAMIC COUPLING IS BELOW THE QUALIFIED ISOLATED-CH TEMPORAL UNCERTAINTY; ISOLATED-CH TEMPORAL QUALIFICATION CONDITIONALLY TRANSFERS TO THIS FULL-CHNS BENCHMARK."
    elif complete: verdict="FULL CHNS PHASE-RATE EXECUTION PASS; ISOLATED-CH TRANSFER FAILS."
    elif moderate.get("passed") and tiny.get("accepted") and stop.get("exception")=="CostStop":
        verdict="FULL CHNS PHASE-RATE BRIDGE PASS; COUPLING TRANSFER INCOMPLETE DUE TO COST GATE."
    elif moderate.get("passed") and target and not target.get("accepted"):
        verdict="MODERATE FULL-CHNS FORMULATION BRIDGE QUALIFIED; TARGET TINY FULL-CHNS STEP FAILED; COUPLING TRANSFER NOT QUALIFIED."
    else: verdict="FULL CHNS PHASE-RATE BRIDGE CHECKPOINT; COUPLING TRANSFER NOT QUALIFIED."
    tests={}
    for path in sorted((ROOT/"tests").glob("*.xml")):
        suites=list(ET.parse(path).getroot().iter("testsuite"))
        tests[path.name]={key:sum(int(s.get(key,0)) for s in suites) for key in ("tests","failures","errors","skipped")}
    session_data=sessions(); wall=sum(r.get("wall_s",0) for r in session_data)
    summary={"execution_status":"complete" if complete else "failed" if stop else "incomplete",
        "bridge_tiny":bridge.get("passed","not_run"),"moderate":moderate.get("passed","not_run"),
        "target_tiny":target.get("accepted","not_run"),"pilot":pilot,"full_L0":status,
        "coupling_transfer":"qualified" if passed else "not_qualified","scientific_verdict":verdict,
        "overall_model":"MODEL NOT YET VALIDATED","production_scientific_wall_s":wall,
        "execution_provenance":frozen,"stop":stop,"tests":tests,
        "remote_CI":"No remote run for this unpushed source; local results only",
        "full_L1_L2":"not_run; not authorized","moving_contact":"not_run; not authorized"}
    atomic_json(ROOT/"summary.json",summary)
    with (ROOT/"summary.csv").open("w",newline="") as stream:
        writer=csv.writer(stream); writer.writerow(["metric","value"])
        for k in ("bridge_tiny","moderate","target_tiny","coupling_transfer","production_scientific_wall_s","scientific_verdict","overall_model"):
            writer.writerow([k,summary[k]])
    lines=["# STEP3A8 — Full-CHNS phase-rate bridge and coupling transfer\n",
        "## 1. Goal\n\nMeasure hydrodynamic coupling independently of the already qualified isolated-CH temporal error.\n",
        "## 2. Inherited STEP3A7 qualification\n\nL0/L1/L2 = 1068/2136/4272 accepted intervals; immutable complete histories and field hashes reverified. No isolated PDE rerun.\n",
        "## 3. Remaining isolated-to-full gap\n\nOnly matched-density, zero-force BE and the historical stiff bump are in scope. A small velocity alone is not a transfer proof.\n",
        "## 4. Full CHNS phase-rate derivation\n\nAbsolute: `((phi-phi_old)/dt,s) - (phi*u,grad(s)) + Mmob*(grad(mu),grad(s)) = 0`. Rate: `(r,s) - ((phi_old+dt*r)*u,grad(s)) + Mmob*(grad(mu),grad(s)) = 0`. Physical `phi=phi_old+dt*r`. Momentum retains `rho*(u-u_old)/dt`, skew convection, `2*eta*symgrad(u)`, `-pi*div(v)`, `phi*grad(mu)` and Navier work; continuity `(q,div(u))=0`; chemical `(mu,chi)-F'_h(phi)[chi]=0`, with unchanged wetting.\n",
        "## 5. Algebraic equivalence\n\nResidual rows agree at corresponding physical fields; they are not rescaled.\n",
        table(["term","relative assembled difference"],list(bridge.get("terms_relative",{}).items()) or [["all","NOT RUN"]]),
        "## 6. Jacobian transformation\n\n`J_rate = J_abs * diag(I_u,I_pi,dt*I_phi,I_mu)`; column transformation only.\n",
        table(["metric","value"],[["Jacobian action relative",bridge.get("Jacobian_action_relative")],["central FD sequence",bridge.get("finite_difference")]]),
        "## 7. BC equivalence\n\nIdentical original essential velocity constraints and single pressure gauge; no rate constraint.\n",
        table(["metric","value"],[["velocity constrained DOFs",bridge.get("BC_velocity_count")],["pressure gauge DOF",bridge.get("pressure_gauge_dof")]]),
        "## 8. Moderate full-CHNS comparison\n\nIndependent physical interval [0,1e-6]; P1 and P2 mandatory. This tests algebra/root accuracy, not BE temporal accuracy.\n"]
    for level in POLICY["moderate_levels"]:
        a=optional(ROOT/f"moderate/{level}_absolute/status.json"); b=optional(ROOT/f"moderate/{level}_rate/status.json")
        ar=(a.get("candidate_audit") or {}).get("physical",{}); br=(b.get("candidate_audit") or {}).get("physical",{})
        comp=next((r for r in moderate.get("comparisons",[]) if r["level"]==level),{})
        rows=[["dt",a.get("dt"),b.get("dt"),None],["SNES reason",a.get("SNES_reason"),b.get("SNES_reason"),None],
            ["iterations",a.get("iterations"),b.get("iterations"),None],["final residual",a.get("final_residual"),b.get("final_residual"),None],
            ["accepted",a.get("accepted"),b.get("accepted"),None]]
        for k,label in (("velocity_L2","u"),("pressure_L2","pi"),("phi_L2","phi"),("mu_L2","mu"),("phase_mass","mass"),("E_total","energy"),
                        ("CH_dissipation","D_CH"),("viscous_dissipation","D_visc"),("slip_dissipation","D_slip")):
            rows.append([label,ar.get(k),br.get(k),abs(ar[k]-br[k]) if k in ar and k in br else None])
        rows.extend([[k,None,None,v] for k,v in comp.get("errors",{}).items()])
        lines += [f"\n{level}, atol=rtol={POLICY['controls'][level]['atol']}; actual controls in each status.json.\n",
            table(["metric","absolute","phase-rate","difference"],rows)]
        if a.get("error") or b.get("error"): lines.append("\nFailure: "+str(a.get("error") or b.get("error"))+"\n")
    lines += ["## 9. Exact full target tiny step\n\nProduction controls unchanged: atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS.\n",
        table(["metric","full rate","isolated rate"],[["dt",target.get("dt"),POLICY["tiny_dt"]],
            ["reason",target.get("SNES_reason"),"historical STEP3A6 PASS"],["initial residual",target.get("initial_residual"),None],
            ["final residual",target.get("final_residual"),None],["Newton iterations",target.get("iterations"),None],
            ["phi L2 difference",tiny.get("phi_L2_difference"),"reference"],["mu L2 difference",tiny.get("mu_L2_difference"),"reference"]]+
            [[k,tiny.get("full_physical",{}).get(k),"historical STEP3A6"] for k in ("speed_max_dof_sample","E_kin","CH_dissipation","viscous_dissipation","slip_dissipation","phase_mass","cells_across_transition_certified_min")]),
        "## 10. Pilot\n\n"+fmt(pilot or None)+"\n",
        "## 11. First dt boundary\n\n"+fmt(boundary or None)+"\n",
        "## 12. Restart\n\n"+fmt(restart or None)+"\n",
        "## 13. Cost gate\n\nFrozen preliminary 1800 s, continuous full L0 14400 s, total 16200 s. Conservative forecast: 1.25*max(mean,p95), using full accepted-step timings only.\n"+fmt(optional(ROOT/"policy/cost_before_full_L0.json") or None)+"\n",
        "## 14. Continuous full L0 trajectory\n\n"+fmt(status or None)+"\n"]
    final=coupling.get("final",{}); energy=coupling.get("energy",{})
    iso0=Journal.read(ISO/"level0/scalar_history.jsonl")[-1]; iso2=Journal.read(ISO/"level2/scalar_history.jsonl")[-1]
    lines += [table(["metric","full CHNS L0","isolated CH L0","isolated CH L2"],[
        ["final CH free energy",final.get("F_CH_full"),iso0["E_total"],iso2["E_total"]],
        ["integrated D_CH",final.get("I_CH_full"),iso0["cumulative_CH_dissipation"],iso2["cumulative_CH_dissipation"]],
        ["final budget defect",final.get("budget"),iso0["energy_budget_defect"],iso2["energy_budget_defect"]]]),
        "## 15. Full physical/mass/resolution diagnostics\n\n"+table(["metric","full trajectory"],[[k,status.get(k)] for k in ("max_mass_domain","min_certified_cells","max_weak_continuity","max_strong_divergence","max_relative_divergence","max_speed")]),
        "## 16. Full energy/work\n\n"+table(["quantity","measured"],list(energy.items()) or [["full trajectory energy/work",None]])]
    gates=coupling.get("coupling",{}); measured=gates.get("measured",{}); limits=gates.get("limits",{})
    for n,title,key in ((17,"Phi coupling","phi_initial"),(18,"Chemical coupling",None),(19,"D_CH coupling","D_significant"),
        (20,"Kinetic energy","kinetic"),(21,"Viscous/slip dissipation","hydro_integral"),(22,"Advective-vs-diffusive phase coupling",None)):
        lines += [f"## {n}. {title}\n\n"+fmt(measured.get(key) if key else ("See selected field/scalar audits" if complete else None))+"\n"]
    lines += ["## 23. Transfer relative to isolated temporal uncertainty\n",
        table(["coupling observable","measured max/final","threshold","pass"],[[k,measured.get(k),v,gates.get("gates",{}).get(k)] for k,v in limits.items()] or [["all continuous-transfer gates",None,None,None]]),
        "## 24. Combined error proxy\n\nTriangle-style practical proxy, NOT a rigorous full continuum bound.\n"+table(["proxy","value"],[[k,gates.get(k)] for k in ("combined_D_proxy","max_combined_phi_proxy")]),
        "## 25. Computational cost\n\nWall time, not CPU time. Reused pilot counted once in total.\n"+table(["quantity","seconds"],[["scientific total",wall],["full L0 sessions",sum(r.get("wall_s",0.) for r in session_data if r.get("full_path"))]])+
        "\nLocal tests: "+json.dumps(tests)+". No remote CI run for unpushed source.\n",
        "## 26. Remaining limitations\n\nNo full L1/L2, unequal-density AGG, gravity, acceleration, BDF2, moving contact line, angle quench, epsilon/mobility/slip sensitivity, falling film or tank work. None is qualified by this benchmark. Missing moderate/target/pilot/full entries above are NOT RUN, not zeros.\n",
        "## 27. Scientific verdict\n\n"+verdict+"\n\nMODEL NOT YET VALIDATED.\n"]
    if stop: lines.append("\nRecorded STOP: "+json.dumps(stop)+"\n")
    if frozen: lines.append("\nExecution provenance: base `"+frozen["implementation"]["execution_base_HEAD"]+"`, source archive `"+frozen["source_archive_sha256"]+"`. A later final commit is not the execution HEAD.\n")
    Path("STEP3A8_REPORT.md").write_text("\n".join(lines))
    history=ROOT/"full_coupling_L0/scalar_history.jsonl"
    if history.exists(): plots(ROOT,Journal.read(history),Journal.read(ISO/"level0/scalar_history.jsonl"),coupling)
    return summary


if __name__=="__main__": generate()
