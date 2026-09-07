"""Standalone scientific PDFs and measured STEP 3A.2 summary; no PDE execution."""
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullLocator

ROOT=Path("validation_results/step3a2")


def read(path):
    return json.loads((ROOT/path).read_text())


def tests(path):
    root=ET.parse(ROOT/path).getroot()
    suites=list(root.iter("testsuite"))
    counts={k:sum(int(s.get(k,0)) for s in suites) for k in ("tests","failures","errors","skipped")}
    counts["passed"]=counts["tests"]-counts["failures"]-counts["errors"]-counts["skipped"]
    counts["green"]=counts["tests"]>0 and counts["failures"]==counts["errors"]==0
    counts["path"]=str(ROOT/path)
    counts["passed_step3a2_cases"]=[case.get("name") for case in root.iter("testcase")
        if "step3a2" in case.get("classname","") and not any(case.find(k) is not None for k in ("failure","error","skipped"))]
    return counts


def history(path):
    with (ROOT/path/"history.csv").open() as stream:
        return list(csv.DictReader(stream))


def main():
    historical=read("p2_resolution/historical_reaudit/actual_coefficients.json")
    prepared=[read(f"equilibrium60/{p}/summary.json") for p in ("prepared_symmetric","prepared_symmetric_fine")]
    metadata=[read(f"equilibrium60/{p}/prepared.json") for p in ("prepared_symmetric","prepared_symmetric_fine")]
    preservation=[read(f"equilibrium60/preservation_{s}/summary.json") for s in ("be","bdf2")]
    convergence=read("equilibrium_perturbation/convergence.json")
    test_results={"pure":tests("tests_pure_final.xml"),"docker":tests("tests_docker_final.xml"),"remote_github_actions_executed":False}
    p=prepared[0]
    pure_cases=test_results["pure"]["passed_step3a2_cases"]
    gates={"P2_extrema_tests":test_results["pure"]["green"] and "test_random_polynomials_never_escape_exact_bounds" in pure_cases,
        "certified_width_tests":test_results["pure"]["green"] and "test_random_gradient_hulls_never_underestimate_width" in pure_cases,**p["gates"],
        "wall_residual_decreases_on_refinement":prepared[1]["prepared"]["wall_residual"]["L2"]<p["prepared"]["wall_residual"]["L2"],
        "BE_preservation":preservation[0]["qualification_status"]=="passed",
        "BDF2_preservation":preservation[1]["qualification_status"]=="passed",
        "perturbation_convergence":convergence["qualification_status"]=="passed",
        "FEniCSx_CI_equivalent":test_results["docker"]["green"]}
    verdict=("DISCRETE EQUILIBRIUM AND BASIC CH TRANSIENT VALIDATED; MOVING-CONTACT STEP 3A TESTS STILL REQUIRED"
             if all(gates.values()) else "MODEL NOT YET VALIDATED")
    comparison=[]
    getters={"phase_mass":lambda r:r["phase_mass"],"apparent_angle_deg":lambda r:r["geometry"]["fits"][1]["theta_deg"],
        "wall_residual_L2":lambda r:r["wall_residual"]["L2"],"wall_residual_max_sampled":lambda r:r["wall_residual"]["max_sampled"],
        "stationary_Riesz_L2":lambda r:r["stationarity"]["L2"],"stationary_Riesz_normalized":lambda r:r["stationarity"]["normalized"],
        "mu_coefficient_spread_Pa":lambda r:r["mu_coefficient_spread"],"D_CH_0_W_per_m":lambda r:r["CH_dissipation"],
        "free_energy_J_per_m":lambda r:r["E_total"],"certified_transition_cells":lambda r:r["resolution"]["cells_across_transition_certified_min"]}
    for metric,get in getters.items():
        comparison.append({"metric":metric,"analytic60_tanh":get(p["analytic60"]),"prepared_equilibrium":get(p["prepared"])})
    result={"scientific_verdict":verdict,"step3a2_pass":all(gates.values()),"gates":gates,
        "historical_interpretation":"geometrically compatible, variationally unprepared",
        "historical_resolution":{"Bernstein_false_positives":historical["bernstein_false_positive_count"],
            "certified_min":historical["exact_p2_certified_normal"]["cells_across_transition_certified_min"],
            "source":"p2_resolution/historical_reaudit/actual_coefficients.json"},
        "comparison_table":comparison,"equilibrium_mesh_study":[{"mesh":r["mesh"],"metadata":m,
            "state":r["prepared"]} for r,m in zip(prepared,metadata)],
        "preservation_table":[{"scheme":r["scheme"],"dt":r["dt"],"t_end":r["t_end"],**r["metrics"]} for r in preservation],
        "perturbation_temporal_table":convergence["table"],"energy_convergence":convergence,
        "read_only_energy_diagnostics":{"initial_rate":read("equilibrium_perturbation/initial_rate_audit.json"),
            "last_BE_intervals":read("equilibrium_perturbation/be_work_audit.json"),"diagnostic_only":True},
        "tests":test_results,"historical_results_modified":False,
        "limitations":["Only matched-density, zero-body-force benchmark at theta_e=60 degrees",
            "Two meshes do not establish a full continuum/epsilon convergence study",
            "No new scientific 90-to-60 spreading, tank, film, water-air or gravity series",
            "Prepared augmented solver/checkpoint currently requires a serial MPI partition",
            "Early failed-preparation summaries lack complete runtime provenance; their source archives and Newton logs are retained",
            "BE work decomposition covers only saved last intervals; continuum trapezoidal budget remains the gate"]}
    (ROOT/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    with (ROOT/"summary.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=comparison[0].keys());writer.writeheader();writer.writerows(comparison)
    with PdfPages(ROOT/"equilibrium_report.pdf") as pdf:
        fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
        axes[0].bar(["analytic tanh","prepared"],[p["analytic60"]["CH_dissipation"],p["prepared"]["CH_dissipation"]])
        axes[0].set(yscale="log",ylabel="initial CH power [W/m]",title="Initial CH impulse removed")
        h=[r["mesh"]["h_max_m"] for r in prepared]
        axes[1].loglog(h,[r["prepared"]["wall_residual"]["L2"] for r in prepared],"o-",label="prepared")
        axes[1].loglog(h,[r["analytic60"]["wall_residual"]["L2"] for r in prepared],"o-",label="analytic FE tanh")
        axes[1].set(xlabel="mesh diameter [m]",ylabel="wall residual L2 [N/sqrt(m)]");axes[1].legend()
        axes[2].plot(h,[r["prepared"]["geometry"]["fits"][1]["theta_deg"] for r in prepared],"o-")
        axes[2].axhline(60.,color="k",ls="--");axes[2].set(xlabel="mesh diameter [m]",ylabel="apparent angle [degree]")
        pdf.savefig(fig);plt.close(fig)
        fig,ax=plt.subplots(figsize=(11,5));ax.axis("off")
        table=ax.table(cellText=[[r["metric"],f'{r["analytic60_tanh"]:.10g}',f'{r["prepared_equilibrium"]:.10g}'] for r in comparison],
            colLabels=["metric","analytic 60-degree tanh","prepared equilibrium"],loc="center",cellLoc="left")
        table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.7)
        ax.set_title("Same inventory, same FE space, same free energy and wetting law")
        pdf.savefig(fig,bbox_inches="tight");plt.close(fig)
    with PdfPages(ROOT/"resolution_report.pdf") as pdf:
        fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        values=[historical["old_centroid_indicator"]["cells_across_transition_min"],
                historical["step3a1_bernstein_sampled"]["cells_across_transition_min"],
                historical["exact_p2_certified_normal"]["cells_across_transition_certified_min"]]
        axes[0].bar(["old centroid","Bernstein / sampled","exact P2 / certified"],values)
        axes[0].axhline(8,color="r",ls="--",label="hard threshold 8");axes[0].legend()
        axes[0].tick_params(axis="x",labelrotation=15);axes[0].set(ylabel="nominal transition cells",title="Historical t=0.5: still FAIL")
        worst=historical["exact_p2_certified_normal"]["worst_cell"]
        triangle=np.array(worst["vertices"])
        axes[1].plot(*np.vstack((triangle,triangle[0])).T,"k-")
        for side,color in (("min","blue"),("max","red")):
            point=worst[side+"_witness"]["physical_point"]
            axes[1].scatter(*point,color=color,label=f'{side}: {worst[side+"_witness"]["value"]:.8f}')
        axes[1].legend();axes[1].set(aspect="equal",xlabel="x [m]",ylabel="z [m]",title="Actual quadratic extrema, cell 8612")
        pdf.savefig(fig);plt.close(fig)
    with PdfPages(ROOT/"energy_report.pdf") as pdf:
        fig,axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
        table=convergence["table"]
        dt=np.array([r["dt"] for r in table])
        axes[0].loglog(dt,[abs(r["final_budget_defect"]) for r in table],"o-")
        axes[0].set(xlabel="dt [s]",ylabel="absolute final budget defect [J/m]")
        axes[1].plot(dt,[r["integrated_D_CH"] for r in table],"o-")
        axes[1].set(xlabel="dt [s]",ylabel="integrated CH power [J/m]")
        axes[2].loglog(dt[:-1],[r["end_phi_difference_to_finest"] for r in table[:-1]],"o-")
        axes[2].set(xlabel="dt [s]",ylabel="end phi L2 difference to finest")
        for ax in axes:
            ticks=dt[:-1] if ax is axes[2] else dt
            ax.set_xticks(ticks,labels=[f"{value:g}" for value in ticks])
            ax.xaxis.set_minor_locator(NullLocator())
        fig.suptitle("Perturbation: decreasing errors, but energy/dissipation qualification FAIL")
        pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        for r in table:
            path=Path(r["run"]).relative_to(ROOT)
            rows=history(path)
            t=np.array([float(row["time"]) for row in rows])
            power=np.array([float(row["CH_dissipation"]) for row in rows])
            delta=np.array([float(row["E_total"]) for row in rows]);delta-=delta[0]
            axes[0].plot(t,power,".-",label=f'dt={r["dt"]:g}')
            axes[1].plot(t,delta,".-",label=f'dt={r["dt"]:g}')
        axes[0].set(xlabel="time [s], including t=0",ylabel="CH power [W/m]")
        axes[1].set(xlabel="time [s]",ylabel="E(t)-E(0) [J/m]")
        for ax in axes:ax.legend()
        pdf.savefig(fig);plt.close(fig)
    print(json.dumps({"verdict":verdict,"gates":gates,"tests":test_results},indent=2))


if __name__=="__main__":
    main()
