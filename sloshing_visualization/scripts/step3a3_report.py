"""Read measured STEP3A3 artifacts; generate a truthful checkpoint report.

No PDE execution and no writes to STEP3/3A1/3A2. Missing nonlinear results are
explicitly NOT RUN, never filled using a linear or low-mode surrogate.
"""
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sloshing.multiphase.ch_timestep_planner import ReducedCHReference
from sloshing.multiphase.benchmarks.spectral_ch import qualified_reference

ROOT = Path("validation_results/step3a3")
OLD = Path("validation_results/step3a2/equilibrium_perturbation")


def read(path):
    return json.loads((ROOT/path).read_text())


def save(path, data):
    (ROOT/path).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def numbers(path):
    with path.open() as stream:
        return {k: np.array([float(r[k]) for r in rows]) for rows in [list(csv.DictReader(stream))] for k in rows[0]}


def test_counts(path):
    root = ET.parse(ROOT/path).getroot()
    suites = list(root.iter("testsuite"))
    result = {k: sum(int(s.get(k, 0)) for s in suites) for k in ("tests", "failures", "errors", "skipped")}
    result["passed"] = result["tests"]-result["failures"]-result["errors"]-result["skipped"]
    result["duration_s"] = sum(float(s.get("time", 0)) for s in suites)
    result["artifact"] = str(path)
    return result


def table(headers, rows):
    def cell(value):
        if value is None:
            return "NOT RUN"
        if isinstance(value, (float, np.floating)):
            return f"{value:.8g}"
        return str(value).replace("|", "\\|")
    return "\n".join(["| "+" | ".join(cell(h) for h in headers)+" |", "| "+" | ".join(["---"]*len(headers))+" |"]+
        ["| "+" | ".join(cell(v) for v in row)+" |" for row in rows])


def main():
    operator = read("linear_operator/status.json")
    spectrum = read("spectrum/spectrum.json")
    ritz_audit = read("spectrum/full_ritz_residual_audit.json")
    if ritz_audit["qualification_status"] != "passed":
        raise ValueError("Actual full Ritz residual audit did not pass")
    for entry in ritz_audit["rows"]:
        row = next(r for r in spectrum["ritz"] if r["mode"] == entry["mode"])
        row["archived_recurrence_estimate"] = row["relative_Ritz_residual"]
        row["relative_Ritz_residual"] = entry["full_operator_relative_residual"]
    rational = read("krylov/rational_balanced/status.json")
    last = read("krylov/rational_balanced/level_m64.json")
    low = read("low_mode/status.json")
    initial_audit = read("linear_operator/initial_audit/status.json")
    planner = read("stiff_bump/planner/status.json")
    plan = read("stiff_bump/planner/timestep_plan.json")
    if planner["reason"] != "STOP: predeclared production cost guard exceeded":
        raise ValueError("This report template is for the measured cost-stop checkpoint, not an arbitrary outcome")
    op, _, _, _, basis, reference = qualified_reference()
    V = basis["V"]
    chemical = op.mass_solve(op.H @ V)
    reduced = ReducedCHReference(basis["projected_L"], V.T @ (op.H @ V),
        op.mobility*chemical.T @ (op.K @ chemical), float(basis["beta"])*np.eye(V.shape[1])[:, 0], V.T @ (op.M0 @ V))
    tests = {"pure": test_counts(Path("tests/pure.xml")), "Docker_CI_equivalent": test_counts(Path("tests/docker.xml")),
        "remote_GitHub_Actions_executed": False,
        "note": "No push or remote workflow dispatch. One incorrect local filename invocation ran zero tests; its XML is retained separately."}
    significant = [r for r in spectrum["ritz"] if r["converged"] and r["significant"]]
    representative = [spectrum["ritz"][0], min(spectrum["ritz"], key=lambda r: abs(np.log(abs(r["real"])/1e7))),
        min(significant, key=lambda r: abs(r["real"])), max(significant, key=lambda r: abs(r["real"]))]
    physical_integral = 1e-6*reduced.integrated_power(1e-4)
    block_table = [[b["block"], b["t_start"], b["t_end"], b["dt"], b["steps"], b["predicted_max_error"],
                    b["active_stiffness"], b["active_stiffness"]*b["dt"]] for b in plan["blocks"]]
    nonlinear = [{"level": n, "status": "NOT RUN: cost STOP", "min_dt": plan["dt0"]/2**n,
        "planned_steps": plan["total_steps"]*2**n, "accepted_steps": 0, "integrated_D_CH": None,
        "DeltaE": None, "final_budget_defect": None, "defect_over_DeltaE": None, "end_phi_difference": None} for n in range(3)]
    old_integrals = [json.loads((OLD/name/"summary.json").read_text())["final"]["cumulative_CH_dissipation"]
        for name in ("be_dt_1e-5", "be_dt_5e-6", "be_dt_2_5e-6")]
    summary = {"scientific_verdict": "MODEL NOT YET VALIDATED", "step3a3_pass": False,
        "stop_reason": planner["reason"], "operator": operator,
        "actual_full_Ritz_residual_audit": ritz_audit, "nonlinear_initial_state_audit": initial_audit,
        "spectrum_summary": {"unrestarted_reference_status": spectrum["qualification_status"],
            "unrestarted_last_error": spectrum["last_reference_comparison"],
            "full_ritz_real_range": [min(r["real"] for r in spectrum["ritz"]), max(r["real"] for r in spectrum["ritz"])],
            "converged_significant_count": len(significant),
            "converged_significant_real_range": [min(r["real"] for r in significant), max(r["real"] for r in significant)],
            "rho_relevant": spectrum["rho_relevant"], "tau_fast_relevant": spectrum["tau_fast_relevant"],
            "max_significant_converged_imaginary": max(abs(r["imag"]) for r in significant),
            "representatives": representative},
        "qualified_reference": {k: rational[k] for k in ("qualification_status", "levels", "dimension", "last_reference_comparison",
            "Krylov_sha256", "reference_sha256", "runtime_s")},
        "reference_energy_identity_max_relative_defect": last["max_projected_energy_identity_relative_defect"],
        "linear_reference_integrated_D2_at_delta_1e_3": physical_integral,
        "low_mode": low, "planner": planner, "planned_block_table": plan["blocks"],
        "linear_plan_power_units_note": "Plan q=psi is normalized by delta; multiply its E2/D2/integrals by delta^2=1e-6 for the unchanged physical bump amplitude",
        "full_nonlinear_levels": nonlinear,
        "full_nonlinear_first_BE_work": read("be_work/stiff_bump_status.json"),
        "old_Richardson_audit_only": {"historical_integrals": old_integrals,
            "extrapolated_last_pair": 2*old_integrals[-1]-old_integrals[-2], "qualification_evidence": False},
        "tests": tests, "historical_results_modified": False,
        "gates": {"operator": True, "relevant_fast_spectrum": True, "qualified_Krylov_reference": True,
            "linear_energy_identity": True, "low_mode": True, "reduced_schedule_accuracy": True,
            "production_cost": False, "full_sparse_schedule_verification": None,
            "full_nonlinear_stiff_bump": None, "nonlinear_energy_dissipation_endpoint_convergence": None}}
    save("summary.json", summary)
    csv_rows = [["operator", "relative mass defect", operator["matrix_checks"]["max_relative_mass_defect"], "PASS"],
        ["operator", "H symmetry", operator["matrix_checks"]["frobenius_symmetry_defects"]["H"], "PASS"],
        ["reference", "max relative L2 self-difference", rational["last_reference_comparison"]["max_relative_L2_difference"], "PASS"],
        ["planner", "level0 predicted hours", plan["cost_estimate"]["level0_hours"], "COST STOP"]]
    csv_rows += [["low_mode", "n="+str(r["steps"])+" final budget/actual change", r["final_budget_over_actual_change"], "PASS"] for r in low["temporal"]]
    csv_rows += [["nonlinear_level"+str(n), "integrated D_CH", "", "NOT RUN"] for n in range(3)]
    with (ROOT/"summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["stage", "metric", "value", "status"]); writer.writerows(csv_rows)

    # Spectrum figure: distinguish converged seeded Ritz data from approximate
    # Ritz values, so a slow unconverged value is never presented as an eigenmode.
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    rates = np.array([-r["real"] for r in spectrum["ritz"]])
    axes[0].loglog(rates, [r["seed_L2_significance"] for r in spectrum["ritz"]], ".", color=".6", label="all m=1024 Ritz (many unconverged)")
    axes[1].loglog(rates, [r["relative_Ritz_residual"] for r in spectrum["ritz"]], ".", color=".6")
    for ax, key in zip(axes, ("seed_L2_significance", "relative_Ritz_residual")):
        ax.loglog([-r["real"] for r in significant], [r[key] for r in significant], "r*", ms=5,
                  label="converged and significant")
        for dt, color in zip((1e-5, 5e-6, 2.5e-6), ("C0", "C1", "C2")):
            ax.axvline(1/dt, color=color, ls=":", label=f"old |lambda| dt=1, dt={dt:g}")
        ax.axvline(.05/plan["dt0"], color="k", ls="--", label="new candidate z=0.05 line")
    axes[0].set(ylabel="seed modal L2 significance", title="Seeded M-Arnoldi Ritz spectrum; reference uses converged rational Krylov")
    axes[1].axhline(1e-6, color="r", ls=":")
    axes[1].set(xlabel="-Re(lambda) [1/s]", ylabel="Ritz residual: recurrence gray, full operator red")
    axes[0].legend(fontsize=7, loc="best")
    fig.savefig(ROOT/"spectrum.pdf"); fig.savefig(ROOT/"spectrum.png", dpi=160); plt.close(fig)

    histories = {name: numbers(OLD/name/"history.csv") for name in ("be_dt_1e-5", "be_dt_5e-6", "be_dt_2_5e-6")}
    t = reference["times"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    axes[0, 0].plot(t, 1e-6*reference["power"], "k-", label="qualified linear D2, delta=1e-3")
    axes[0, 1].plot(t, [1e-6*reduced.integrated_power(x) for x in t], "k-", label="integral of linear D2")
    axes[0, 1].plot(t, 1e-6*(reference["energy"][0]-reference["energy"]), "--", label="linear E2(0)-E2(t)")
    for name, h in histories.items():
        label = "historical "+name.replace("be_dt_", "dt=")
        axes[0, 0].plot(h["time"], h["CH_dissipation"], ".", label=label+" snapshots")
        axes[0, 1].plot(h["time"], h["cumulative_CH_dissipation"], ".--", label=label)
        axes[1, 0].plot(h["time"], h["energy_budget_defect"], ".--", label=label)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 0]):
        ax.set_xscale("symlog", linthresh=t[1]); ax.set_xlabel("time [s]; symlog retains t=0")
        ax.set_xticks([0., 1e-12, 1e-10, 1e-8, 1e-6, 1e-4],
                      labels=["0", "1e-12", "1e-10", "1e-8", "1e-6", "1e-4"])
        ax.legend(fontsize=7)
    axes[0, 0].set(yscale="log", ylabel="CH power [W/m]")
    axes[0, 1].set(ylabel="cumulative physical power / energy loss [J/m]")
    axes[1, 0].set(ylabel="historical budget defect [J/m]")
    axes[1, 1].axis("off")
    axes[1, 1].text(.04, .85, "NEW FULL NONLINEAR LEVELS: NOT RUN\n\nCost STOP: 4.65 h level0; 32.54 h series.\nNo new nonlinear curve is inferred from\nlinear reference or isolated low mode.\n\nMODEL NOT YET VALIDATED", va="top", fontsize=11)
    fig.suptitle("Stiff startup: linear diagnosis and unchanged historical observations")
    fig.savefig(ROOT/"stiff_startup_energy.pdf"); plt.close(fig)

    first_folder = ROOT/"low_mode/be_n20_maxdelta_0.0001"
    audits = [json.loads(line) for line in (first_folder/"be_work.jsonl").read_text().splitlines()]
    save("be_work/low_mode_early_intervals.json", {"scope": "isolated CH low mode, NOT the stiff bump", "intervals": audits[:20]})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for key, sign, label in (("weak_work_defect", 1, "weak-work defect"), ("B_BE", -1, "-B_BE (numerical)"),
                            ("trapezoid_gap", 1, "dt/2 (D_old-D_new)"), ("continuum_local_defect", 1, "continuum local defect")):
        axes[0].plot([r["interval"] for r in audits], [sign*r[key] for r in audits], ".-", label=label)
    axes[0].set(yscale="symlog", xlabel="BE interval", ylabel="signed work [J/m]", title="Measured isolated CH low mode")
    axes[0].set_yscale("symlog", linthresh=1e-18)
    axes[0].legend(fontsize=8, loc="upper left", bbox_to_anchor=(0., -.15), ncol=2)
    axes[1].axis("off"); a = audits[0]
    axes[1].text(.02, .95, "ISOLATED LOW MODE, interval 1\n\n"+
        "\n".join(f"{k}: {a[k]:.8e} J/m" for k in ("weak_work_defect", "B_BE", "trapezoid_gap", "continuum_local_defect"))+
        "\n\ncontinuum = weak - B_BE + trapezoid gap\nB_BE is NOT physical heat.\n\nNew stiff-bump intervals: NOT RUN.", va="top", fontsize=10)
    fig.savefig(ROOT/"early_be_work_decomposition.pdf"); plt.close(fig)

    checks = operator["matrix_checks"]
    finite = operator["smooth_bump_finite_differences"]
    low_rows = [[r["steps"], r["dt"], r["measured_modal_amplitude"], r["linear_BE_amplitude"], r["error_to_exponential"],
        r["integrated_D_CH"], r["DeltaE"], r["final_budget_defect"], r["final_budget_over_actual_change"]] for r in low["temporal"]]
    report = f"""# STEP 3A.3 — CH spectral validation checkpoint

**MODEL NOT YET VALIDATED.** Operator, rational exponential reference and
isolated low-mode checks PASS. Full nonlinear stiff-bump qualification is
**NOT RUN**: the predeclared production cost gate stopped execution.

## 1. Goal

Resolve the CH-dominated fast spectrum without changing the model, initial
equilibrium, P2 space, mesh, quadrature, mobility or saved delta=1e-3 bump.
The sequence and numerical thresholds were fixed in [STEP3A3_DESIGN.md](STEP3A3_DESIGN.md)
before PDE runs. No first interval or t=0 power is omitted. Historical reports
and step3/step3a1/step3a2 artifacts are unchanged.

## 2. STEP3A.2 unresolved temporal problem

Initial excess energy was 5.542888467657825e-8 J/m and D_CH(0) was
0.11400313905309394 W/m. Thus tau_E=4.862048987156724e-7 s; old dt/tau_E
were {', '.join(f'{x:.8g}' for x in operator['old_dt_over_tau_E'])}.
The old nonzero transient FAIL is not rewritten. Its roundoff-scale initial
semidiscrete and last-interval BE identities motivate a time-resolution
investigation, not a free-energy or weak-form change.

## 3. Linearized CH derivation

M q_dot + mobility K r = 0, M r = H q, hence

`L = -mobility M^(-1) K M^(-1) H`.

Only sparse matrices and solves are used at full size. H is obtained with UFL
as the second derivative of the existing discrete scalar free energy.
This analysis does not replace the production CHNS solver.

## 4. Discrete matrices M/K/H

There are {operator['phase_dofs']} phase DOFs, degree-12 quadrature, P2 and the
unchanged 96x48 crossed mesh. H includes bulk
`lambda/epsilon*(3*phi_e^2-1) q chi`, the gradient term and the exact wall term.
Current `f_wall=-sigma*cos(theta)*(3phi-phi^3)/4`, so
`f_wall''=3*sigma*cos(theta)*phi/2`; no wall Hessian approximation is used.

{table(['matrix', 'Frobenius symmetry defect'], [[k, v] for k, v in checks['frobenius_symmetry_defects'].items()])}

## 5. Mass-zero projection

`Pq=q-1_h*(m^T q)/(m^T 1_h)`, with assembled `m_i=integral N_i`.
Projection is permitted only inside analysis; there is no production mass fix.
Twenty deterministic vectors give max relative mass defect
{checks['max_relative_mass_defect']:.10g} and relative `PLPq-Lq` defect
{checks['max_projection_relative_defect']:.10g}. Both pass 1e-11.

## 6. Finite-difference operator validation

An independent nonlinear weak chemical assembly produces the semidiscrete rate,
without time advance. Central differences remove the O(delta) nonlinear term.
The one-sided limit is reported too, not confused with central accuracy.

{table(['delta', 'Hessian M-dual relative defect', 'central rate relative L2', 'one-sided rate relative L2'], [[r['delta'], r['hessian_relative_M_dual_error'], r['central_rate_relative_L2_error'], r['one_sided_rate_relative_L2_error']] for r in finite])}

Three additional deterministic random FE directions also PASS. The final
central-rate plateau is roundoff; the mandatory 1e-6 operator gate passes.

## 7. CH-dominated early-time justification

Read-only historical early samples give max `(D_visc+D_slip)/D_CH` =
{operator['CH_dominated_justification']['max_hydrodynamic_power_over_CH']:.10g}
and max `E_kin/initial_excess` =
{operator['CH_dominated_justification']['max_kinetic_over_initial_excess']:.10g}.
The spectral planner analyzes the CH-dominated stiff subsystem, while final
qualification remains a FULL nonlinear CHNS experiment, not performed here.

## 8. Relevant Ritz spectrum

The own M-inner-product, twice-reorthogonalized Arnoldi was seeded with the
historical psi coefficients. Dimensions 32/48/64 were increased through 1024.
All {len(significant)} converged significant Ritz values are real and negative;
their range is {summary['spectrum_summary']['converged_significant_real_range']} s^-1.
Not every approximate Ritz value is a converged eigenmode. The full spectrum
is not claimed resolved from this subset.

Initial saved residuals were terminal Arnoldi recurrence estimates. An
independent actual full-operator audit of all 33 converged significant vectors
passed, max residual {ritz_audit['max_full_operator_relative_residual']:.10g}.
It replaces sub-roundoff estimates in the red plot markers and converged table
rows, without rewriting the saved primary spectrum. Unconverged table rows
still show labelled recurrence estimates. The current API uses full residuals.

{table(['Ritz mode', 'Re(lambda) [1/s]', 'Im', 'tau [s]', 'relative residual', 'old z(2.5e-6)', 'new z(dt0)', 'converged'], [[r['mode'], r['real'], r['imag'], r['tau_s'], r['relative_Ritz_residual'], r['old_z_finest'], abs(r['real'])*plan['dt0'], r['converged']] for r in representative])}

Complete records: `validation_results/step3a3/spectrum/ritz.csv`.
Figure: [spectrum.pdf](validation_results/step3a3/spectrum.pdf),
[spectrum.png](validation_results/step3a3/spectrum.png).

## 9. Krylov self-convergence

The single polynomial reference FAILED: m=768/1024 max relative L2 difference
{spectrum['last_reference_comparison']['max_relative_L2_difference']:.8g}.
This is not a detected physical instability: orthogonality and the Arnoldi
relation remained near machine precision. `rho*t_end` is about 1.07e6.

The installed SLEPc 3.24 MFN alternative was investigated without installing
packages. A tiny nonidentity-M test rejected direct BV/M normalization; the
correct Euclidean-internal adapter passed its dense test. Its full-horizon
ncv=32 run reached t=7.711254e-5 at 2736.52 s and was interrupted at the
one-hour resource guard before completing the horizon. It is NOT qualified.
See the primary [SLEPc MFN documentation](https://slepc.upv.es/release/documentation/manual/mfn.html)
and [SciPy matrix exponential documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html)
for the underlying tooling; numerical evidence here comes from saved runs.

A separately versioned rational M-Krylov reference passed levels 32/48/64.
Fixed positive resolvent times are geometrically spaced from 1/rho to 1e-4 s.
The sparse BE resolvent only generates a basis; reference evolution is
`V exp(t V^T M L V) e1 * ||q0||_M`, NOT a BE time history. Balanced auxiliary
scaling and pinned MUMPS change only linear algebra. No model or gate changed.

{table(['48/64 comparison', 'maximum relative difference', 'threshold'], [[k, v, 1e-6] for k, v in rational['last_reference_comparison'].items() if k != 'relative_L2_by_time'])}

M-orthogonality defect is {last['diagnostics']['orthogonality_defect']:.10g}.
The rational reference's individual Ritz vectors are not asserted converged;
its exponential action is self-converged in full FEM L2/E2/D2. Fast relevant
rates for the planner come from the independently converged primary Ritz set.

## 10. Semidiscrete linear energy identity

`E2=q^T H q/2`, `D2=mobility*(M^-1 Hq)^T K(M^-1 Hq)`.
The full operator's max relative `dE2/dt+D2` defect is
{checks['max_relative_energy_identity_defect']:.10g}; the qualified reference's
projected identity defect is {last['max_projected_energy_identity_relative_defect']:.10g}.
Both pass the 1e-10 identity gate. Integrated reference power is evaluated
independently from its quadratic modal integrand, never replaced with -DeltaE.

## 11. Initial spectral timescales

At the unchanged delta=1e-3, linear E2(0)={operator['initial_linear_E2_at_delta']:.16g}
J/m, D2(0)={operator['initial_linear_D2_at_delta']:.16g} W/m and
tau_E_linear={operator['tau_E_linear']:.16g} s. This agrees with the historical
nonlinear timescale. The fastest converged significantly seeded rate is
rho={spectrum['rho_relevant']:.16g} s^-1, tau={spectrum['tau_fast_relevant']:.16g} s.
Its weight is small, but exceeds the predeclared significance criterion; this
is not an unseeded all-grid spectral-radius estimate. Energy timescale and the
fastest weakly excited timescale are different quantities.

## 12. Low-mode temporal benchmark

An independently refined full-operator mode has lambda={low['mode']['real']:.16g}
s^-1, tau={low['mode']['tau']:.16g} s, residual
{low['mode']['relative_full_operator_residual']:.10g}. Its M norm is 1;
max initial |delta phi|=1e-4. This different input is explicitly the isolated
low-mode benchmark, NOT a replacement for the fixed stiff bump.

{table(['steps', 'dt', 'measured amplitude', 'linear BE amplitude', 'error to exp(-1)', 'integral D_CH', 'DeltaE', 'budget defect', 'defect/|DeltaE|'], low_rows)}

Observed orders: {low['observed_modal_orders']}. Amplitude-halving relative BE
errors: {low['amplitude_BE_error_sequence']}. Max mass/domain error over all
five low-mode runs is {max(r['max_mass_domain_error'] for r in low['temporal']+low['amplitude_checks']):.10g};
minimum certified resolution is {min(r['minimum_certified_cells'] for r in low['temporal']+low['amplitude_checks']):.10g}.
All isolated low-mode gates PASS, including integrated D and continuum budget
convergence. This does not qualify full CHNS on the stiff bump.

## 13. Spectral timestep planner

z_target=0.05 gives initial candidate dt0={plan['candidate_dt0']:.16g} s.
The candidate is proposed entirely from the saved linear input; blocks grow
by at most two, with a 1e-4 initial-D2 content cutoff used only for proposal.
The first proposal failed strict global linear error gates. Uniform halving
passed: dt0={plan['dt0']:.16g} s, first-block z={plan['dt0_over_tau_fast']},
single-step BE modal amplification error={plan['max_BE_modal_amplification_error_first_block']:.10g}.
Both candidates, including the rejection, are retained. No nonlinear observed
field can modify a frozen schedule, and an unauthorized plan is rejected by
the execution-schedule API.

## 14. Planned schedule

This is an immutable **reduced-linear-verified candidate, NOT an authorized
production plan**. Total 1796 intervals to 1e-4 s; half/quarter plans would use
3592/7184 intervals, identical boundaries. `active z` uses content still above
the proposal cutoff; the complete JSON also records rho_initial*dt for every
block (including already-decayed modes).

{table(['block', 't_start', 't_end', 'dt', 'steps', 'predicted max L2 error', 'active stiffness [1/s]', 'active z'], block_table)}

File: [timestep_plan.json](validation_results/step3a3/stiff_bump/planner/timestep_plan.json).
File SHA256: `{planner['timestep_plan_sha256']}`.
Semantic plan SHA256: `{plan['plan_sha256']}`.

## 15. Linear BE schedule verification

{table(['predicted reduced error', 'value', 'threshold'], [[k, plan['predicted_errors'][k], limit] for k, limit in [('endpoint_relative_L2_error', 1e-3), ('max_checkpoint_relative_L2_error', 1e-3), ('max_quadratic_energy_relative_error', 1e-3), ('integrated_D2_relative_error', 1e-2)]])}

Reference power integral at delta=1e-3 is {physical_integral:.16g} J/m;
scheduled reduced BE predicts {1e-6*plan['predicted_errors']['linear_BE_integrated_D2']:.16g} J/m.
The plan's stored unscaled q=psi powers/integrals must be multiplied by delta².
Sparse full-size BE verification is implemented and tested against a dense
tiny tangent system, but its scientific schedule run was **NOT RUN after cost STOP**.
Reduced verification alone does not authorize production nonlinear execution.

## 16. Full nonlinear stiff-bump level 0

NOT RUN. Historical conservative seconds/accepted interval =
{plan['cost_estimate']['seconds_per_step']:.10g}, giving
{plan['cost_estimate']['level0_hours']:.10g} hours for level0 and
{plan['cost_estimate']['level0_half_quarter_hours']:.10g} hours for all three
levels. Predeclared limits are 3 and 10 hours. They were not increased after
seeing the candidate. The measured blocker is computational cost, not a newly
observed weak-form defect. A separate zero-time audit reloaded exactly the
frozen field and re-ran only the existing consistent weak mu initialization:
D_CH(0)={initial_audit['initial']['CH_dissipation']:.16g} W/m,
excess energy={initial_audit['initial_excess_energy']:.16g} J/m. The actual
mixed-state phi hash matches the historical initial hash exactly. This audit
has ZERO accepted time intervals and does not bypass the cost STOP.

## 17. Half-plan refinement

NOT RUN: level0 not authorized. The deterministic rule and 3592-interval
schedule exist; no measured full nonlinear refinement claim is made.

## 18. Optional quarter-plan refinement

NOT RUN: level0 not authorized. The deterministic rule and 7184-interval
schedule exist, but cost blocks the series and endpoint/integral order evidence.

## 19. Integrated dissipation convergence

{table(['level', 'min dt (planned)', 'planned steps', 'accepted', 'integral D_CH', 'DeltaE', 'final defect', 'defect/|DeltaE|', 'end phi difference'], [[r['level'], r['min_dt'], r['planned_steps'], r['accepted_steps'], r['integrated_D_CH'], r['DeltaE'], r['final_budget_defect'], r['defect_over_DeltaE'], r['end_phi_difference']] for r in nonlinear])}

The full-CHNS successive D_CH differences and final <=5% gate are unmeasured.
The low-mode PASS and linear integral are not substituted for these gates.

## 20. Energy closure

Full nonlinear final budget, observed budget order and endpoint convergence
remain unmeasured. For audit only, the historical last-pair Richardson hint is
{2*old_integrals[-1]-old_integrals[-2]:.10g} J/m; it is NOT qualification evidence.
The new linear reference is consistent with the hypothesis of temporal error,
but cannot confirm the nonlinear continuum limit by itself.
Figure: [stiff_startup_energy.pdf](validation_results/step3a3/stiff_startup_energy.pdf)
shows the qualified linear reference and labelled historical data only;
t=0 is visible on a symlog axis, not removed.

## 21. Early BE work decomposition

Physical continuum dissipation is **D_CH + D_visc + D_slip**.
**B_BE is a signed numerical/discrete BE work remainder, diagnostic only**.
It can contain negative quartic/cubic pieces and is never clipped or added to
the physical power curve. The identity is
`continuum_local = weak_work_defect - B_BE + dt/2*(D_old-D_new)`.

Measured isolated low-mode first interval (20-step run):

{table(['quantity', 'value'], [[k, a[k]] for k in ('Delta_E', 'D_old', 'D_new', 'trapezoid_integral', 'BE_endpoint_integral', 'weak_work_defect', 'B_BE', 'trapezoid_gap', 'continuum_local_defect', 'decomposition_roundoff')])}

Thus even this resolved low mode has a visible signed BE/trapezoid split while
weak work is near roundoff. Every low-mode interval was audited. A pure stiff
z=100 regression also has D_old/D_new=10201 with a large trapezoid error and
roundoff BE work. New full stiff-bump first-20 audits are **NOT RUN**.
Figure: [early_be_work_decomposition.pdf](validation_results/step3a3/early_be_work_decomposition.pdf)
is labelled isolated CH, never presented as full stiff-bump evidence.

## 22. Linear-vs-nonlinear comparison

Available: t=0 nonlinear/linear timescale agreement and measured low-mode BE
decay versus its linear prediction. Unavailable: common-early-time full-CHNS
bump L2/H1/excess-energy/power comparison. No new nonlinear snapshot exists.
Low mode answers the isolated modal time-integrator question; stiff bump would
answer practical resolution of the broad fast spectrum. They remain separate.

## 23. Remaining limitations

The cost STOP requires a separately agreed compute resource or implementation
cost reduction before full sparse verification and full nonlinear series. No
claim is made that the candidate is runtime-optimal. Full CHNS production
schedule execution and its first-20 work audit remain unperformed, and no
small-dt production conditioning conclusion can be inferred without that run.

No spreading, theta sweep, epsilon/mobility/slip study, film, tank, unmatched
density, gravity or water-air validation was attempted. There is no continuum
or full model qualification here.

Tests: {tests['pure']['passed']} pure passed, {tests['pure']['skipped']} skipped
(FEniCSx/PETSc/SLEPc unavailable on host; 2 expensive tests deselected by the
existing default configuration); {tests['Docker_CI_equivalent']['passed']}
Docker tiny checks passed. Existing negative-resolution smoke emits its
expected warning. `.github/workflows/step3.yml` includes new tiny tests,
never the production 96x48 series. GitHub-hosted Actions were NOT run: no push
or remote workflow dispatch. XML evidence is in `validation_results/step3a3/tests/`.

Provenance: base git `{operator['provenance']['git_commit']}`, pinned image
`{operator['provenance']['docker_image_digest']}`. Each completed stage stores
source hash/archive, equilibrium/mesh/input fingerprints; plan stores policy
and Krylov hashes. Early interrupted reference prototypes retain source
archives and explicit terminal status but not complete final provenance.
The qualified rational basis/reference, sparse matrices and exact recovered
input arrays are retained. The failed 1024-vector basis is local-only, with
all numerical summaries and hashes retained in git.

Prepared equilibrium fingerprint: `{operator['provenance']['equilibrium_fingerprint']}`.
Psi coefficient SHA256: `{operator['provenance']['perturbation_coefficient_sha256']}`.
Initial phi SHA256: `{operator['provenance']['initial_phi_sha256']}`.
Both initial hashes match historical STEP3A.2 exactly. The historical bump
had only hashes, not a standalone saved psi array; recovery replayed the
unchanged initializer once, verified BOTH hashes and froze the identical bytes.

## 24. Scientific verdict

**MODEL NOT YET VALIDATED.**

PASS: sparse operator/finite differences/conservation/symmetry, relevant fast
Ritz evidence, rational exponential reference, linear energy identity,
isolated low-mode temporal convergence and reduced-linear candidate accuracy.

STOP: predeclared production cost limit. NOT RUN: full sparse candidate
verification and full nonlinear stiff-bump plan/half/quarter convergence.
Therefore **BASIC NONZERO CH TRANSIENT TEMPORALLY QUALIFIED** is not claimed.
Historical STEP3A.2 FAIL remains unchanged. No physics or thresholds were
relaxed, and no numerical BE remainder was relabelled as physical heat.
"""
    Path("STEP3A3_REPORT.md").write_text(report)
    save("report_provenance.json", {"generator": "scripts/step3a3_report.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "report_sha256": hashlib.sha256(Path("STEP3A3_REPORT.md").read_bytes()).hexdigest(),
        "data_scope": "read-only measured STEP3A3 plus explicitly labelled historical STEP3A2",
        "scientific_verdict": summary["scientific_verdict"]})
    print(json.dumps({"verdict": summary["scientific_verdict"], "tests": tests, "cost_stop": plan["cost_estimate"]}, indent=2))


if __name__ == "__main__":
    main()
