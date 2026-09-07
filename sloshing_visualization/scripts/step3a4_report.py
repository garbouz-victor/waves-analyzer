"""Measured STEP3A4 checkpoint report. No solves or writes to historical iterations."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sloshing.multiphase.ch_validation_policy import LINEAR_LIMITS, authorizations
from sloshing.multiphase.provenance import source_hash

ROOT = Path("validation_results/step3a4")


def read(path):
    return json.loads((ROOT/path).read_text())


def save(path, data):
    target = ROOT/path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def save_not_run(path, data):
    target = ROOT/path
    if target.exists() and read(path).get("status") != "not_run":
        raise ValueError("Refusing to replace execution evidence with a NOT RUN placeholder")
    save(path, data)


def counts(path):
    suites = list(ET.parse(ROOT/path).getroot().iter("testsuite"))
    record = {k: sum(int(s.get(k, 0)) for s in suites) for k in ("tests", "failures", "errors", "skipped")}
    record["passed"] = record["tests"]-record["failures"]-record["errors"]-record["skipped"]
    record["duration_s"] = sum(float(s.get("time", 0)) for s in suites)
    record["artifact"] = path
    return record


def table(headers, rows):
    def cell(v):
        if v is None:
            return "N/A — NOT RUN"
        return f"{v:.9g}" if isinstance(v, (float, np.floating)) else str(v).replace("|", "\\|")
    return "\n".join(["| "+" | ".join(headers)+" |", "| "+" | ".join(["---"]*len(headers))+" |"]+
        ["| "+" | ".join(cell(v) for v in row)+" |" for row in rows])


def main():
    baseline = read("linear_sparse/baseline_plan/status.json")
    optimized = read("linear_sparse/optimized_plan/status.json")
    optimization = read("schedule_optimization/optimization.json")
    plan = read("schedule_optimization/optimized_plan.json")
    pilot = read("performance/isolated_pilot/status.json")
    input_provenance = read("performance/isolated_pilot/input_provenance.json")
    initial = read("performance/isolated_pilot/initial.json")
    audit = read("performance/failed_interval_audit/status.json")
    selections = read("performance/isolated_pilot/run_plan.json")["diagnostic_selection"]
    if (baseline["qualification_status"] != "passed" or optimized["qualification_status"] != "passed" or
            pilot["status"] != "stopped" or pilot["accepted_steps"] != 0):
        raise ValueError("This checkpoint report requires the measured baseline/optimized PASS and first-interval STOP")
    timing_paths = {"initial_profile": "performance/initial_profile/timing.json",
        "baseline_sparse": "linear_sparse/baseline_plan/timing.json",
        "optimizer_and_repeat": "schedule_optimization/timing.json",
        "optimized_sparse": "linear_sparse/optimized_plan/timing.json",
        "isolated_failed_pilot": "performance/isolated_pilot/timing.json",
        "read_only_failure_audit": "performance/failed_interval_audit/timing.json"}
    timings = {k: read(v) for k, v in timing_paths.items()}
    archive_checks = []
    for value in timing_paths.values():
        folder = Path(value).parent
        metadata = read(folder/"provenance.json")
        path = ROOT/folder/"multiphase_source.tar.gz"
        digest = hashlib.sha256()
        with tarfile.open(path) as archive:
            for member in sorted(archive.getmembers(), key=lambda m: m.name):
                digest.update(member.name.encode())
                digest.update(archive.extractfile(member).read())
        valid = (digest.hexdigest() == metadata["multiphase_source_sha256"] and
                 hashlib.sha256(path.read_bytes()).hexdigest() == metadata["source_archive_sha256"])
        if not valid:
            raise ValueError("Stage source archive does not match provenance")
        archive_checks.append({"stage": str(folder), "source_and_archive_hashes_match": valid})
    if (plan["optimizer_result_sha256"] != optimization["plan_sha256"] or
            any(plan[k] != optimization[k] for k in ("blocks", "predicted_errors", "input_fingerprints", "policy_sha256"))):
        raise ValueError("Frozen execution plan differs from reproducible optimizer core")
    reproducibility = {"independent_optimizer_repeat_matched": plan["optimizer_reproducible"],
        "optimizer_core_semantic_sha256": optimization["plan_sha256"],
        "frozen_execution_artifact_sha256": plan["plan_sha256"],
        "semantics": "Core hash covers deterministic optimizer inputs, policies, history, blocks/errors. Execution hash additionally covers archival provenance; gzip timestamps are not optimizer input.",
        "core_and_executed_schedule_match": True}
    save("schedule_optimization/reproducibility.json", reproducibility)
    old_s = 9.31865964329918
    old_series_h = 7*1796*old_s/3600
    wall = sum(t["wall_s"] for t in timings.values())
    cpu = sum(t["cpu_s"] for t in timings.values())
    analysis_wall = sum(t["wall_s"] for t in timings.values() if t["analysis_budget_charged"])
    provenance = {"report_generation_git_SHA": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "multiphase_source_sha256": source_hash(),
        "report_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "environment": "host report generation; scientific Docker/source provenance is per stage",
        "policy": read("policy/policy.json"),
        "stage_provenance_paths": [str(Path(v).parent/"provenance.json") for v in timing_paths.values()],
        "archive_hash_checks": archive_checks}
    performance = {"historical_full_CHNS_seconds_per_accepted_interval": old_s,
        "historical_full_CHNS_three_level_estimate_wall_hours": old_series_h,
        "measured_wall_s": wall, "measured_process_CPU_s": cpu,
        "analysis_budget_used_wall_s": analysis_wall, "stages": timings,
        "isolated_seconds_per_accepted_step": None,
        "measured_validation_CPU_savings": None,
        "CPU_savings_note": "No accepted isolated step or qualified series. Old 32.54 h is a wall estimate, not measured CPU. No completed-validation speedup claimed.",
        "provenance": provenance}
    save("performance/timing_breakdown.json", performance)
    tests = {"pure": counts("pure_tests.xml"), "pinned_Docker_workflow": counts("docker_tests_completed.xml"),
        "github_remote": "NOT RUN: new commit intentionally not pushed; workflow updated and run locally",
        "initial_Docker_invocation": "docker_tests.xml: no tests ran because host shell did not expand container-relative globs; corrected shell invocation only, not a scientific retry"}
    levels = []
    for level in range(3):
        record = {"status": "not_run", "qualification_status": "not_qualified", "level": level,
            "reason": "STOP: isolated cost pilot rejected first SNES interval; no retry",
            "planned_steps": plan["total_steps"]*2**level, "accepted_steps": 0,
            "min_dt": plan["dt0"]/2**level, "t_end": 1e-4,
            "integrated_D_CH": None, "DeltaF": None, "budget_defect": None,
            "endpoint_L2_difference": None, "optimized_plan_sha256": plan["plan_sha256"]}
        levels.append(record)
        save_not_run(f"isolated_stiff_bump/level{level}/status.json", record)
    probe = {"status": "not_run", "execution_authorized": False,
        "reason": "Isolated nonlinear temporal convergence prerequisite NOT obtained",
        "reference_derived_checkpoint_times": selections["probe_times"],
        "reference_derived_checkpoint_indices": selections["probe_indices"],
        "target_times": selections["probe_target_times"],
        "times_semantics": "predeclared selectors evaluated before pilot; not an authorized/frozen CHNS execution plan",
        "prefix_steps_if_later_authorized": selections["probe_indices"][-1],
        "historical_cost_prefix_estimate_wall_s": selections["probe_indices"][-1]*old_s,
        "optimized_plan_sha256": plan["plan_sha256"], "comparison_metrics": None}
    save_not_run("coupling_probe/planning/status.json", probe)
    for name in ("run_early", "run_mid", "run_late"):
        save_not_run(f"coupling_probe/{name}/status.json", {"status": "not_run", "reason": probe["reason"]})
    result = {"scientific_verdict": "MODEL NOT YET VALIDATED",
        "step3a4_verdict": "LINEAR SCHEDULE VERIFIED; NONLINEAR ISOLATED CH STOPPED AT FIRST SNES INTERVAL; NO FULL-CHNS TRANSFER QUALIFICATION",
        "authorizations": authorizations(reduced=True, sparse=True),
        "baseline_sparse": baseline, "optimized_sparse": optimized,
        "optimized_steps": plan["total_steps"], "percent_shorter_than_1796": 100*(1-plan["total_steps"]/1796),
        "optimized_plan_sha256": plan["plan_sha256"], "optimizer_reproducible": plan["optimizer_reproducible"],
        "reproducibility": reproducibility,
        "initial": initial, "pilot": pilot, "rejected_iterate_audit": audit,
        "isolated_levels": levels, "coupling_probe": probe,
        "last_relative_D_CH_difference": None, "endpoint_L2_differences": None,
        "tests": tests, "timing_summary": {"wall_s": wall, "process_CPU_s": cpu, "analysis_wall_s": analysis_wall},
        "provenance": provenance}
    save("summary.json", result)
    with (ROOT/"summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["stage", "metric", "value", "status"])
        for stage, data in (("baseline_sparse", baseline), ("optimized_sparse", optimized)):
            for key, value in data["errors"].items():
                writer.writerow([stage, key, value, "PASS"])
        for key in ("SNES_reason", "SNES_iterations", "SNES_residual", "accepted_steps"):
            writer.writerow(["isolated_pilot", key, pilot[key], "STOP"])
        for level in levels:
            for key in ("integrated_D_CH", "DeltaF", "budget_defect", "endpoint_L2_difference"):
                writer.writerow([f"isolated_level{level['level']}", key, "", "NOT RUN"])
        writer.writerow(["overall", "scientific_verdict", result["scientific_verdict"], "NOT QUALIFIED"])
    original_errors = optimization["optimization_history"][0]["errors"]
    plans = [("original 898", 898, original_errors, "NOT RUN: reduced FAIL"),
        ("uniform 1796", 1796, baseline["reduced_errors"], "PASS"),
        ("optimized 1068", plan["total_steps"], plan["predicted_errors"], "PASS")]
    plan_table = table(["plan", "steps", "endpoint L2", "max L2", "energy error", "D2 error", "full-sparse pass"],
        [[label, n]+[errors[k] for k in LINEAR_LIMITS]+[status] for label, n, errors, status in plans])
    sparse_table = table(["metric", "reduced optimized", "full sparse optimized", "full minus reduced"],
        [[key, plan["predicted_errors"][key], optimized["errors"][key],
          optimized["errors"][key]-plan["predicted_errors"][key]] for key in LINEAR_LIMITS]+
        [["max mass/domain", "not independently assessed here", optimized["errors"]["max_mass_domain"], "N/A"]])
    blocks_table = table(["block (0-based)", "t_start", "t_end", "dt", "steps"],
        [[b[k] for k in ("block", "t_start", "t_end", "dt", "steps")] for b in plan["blocks"]])
    history_table = table(["candidate", "action", "block", "steps", "endpoint", "max L2", "energy", "D2", "PASS"],
        [[r["candidate"], r["action"], r["block"] if r["block"] is not None else "—", r["steps"]]+
         [r["errors"][k] for k in LINEAR_LIMITS]+[r["passed"]] for r in optimization["optimization_history"]])
    cost_table = table(["stage", "sec/accepted step", "steps", "total wall runtime [s]"],
        [["historical full CHNS (projected 3-level series)", old_s, 7*1796, 7*1796*old_s],
         ["full sparse linear baseline", timings["baseline_sparse"]["wall_s"]/1796, 1796, timings["baseline_sparse"]["wall_s"]],
         ["full sparse linear optimized", timings["optimized_sparse"]["wall_s"]/1068, 1068, timings["optimized_sparse"]["wall_s"]],
         ["isolated failed pilot", "undefined: 0 accepted", "0 accepted / 1 attempted", timings["isolated_failed_pilot"]["wall_s"]]]+
        [[f"isolated CH level{j}", None, f"{1068*2**j} planned / 0 accepted", None] for j in range(3)]+
        [["full CHNS coupling probe", None, "745 prospective / 0 executed", None]])
    series_table = table(["level", "min dt (planned)", "steps planned/accepted", "integral D_CH", "DeltaF", "final defect", "defect/DeltaF", "endpoint difference"],
        [[r["level"], r["min_dt"], f"{r['planned_steps']}/0", None, None, None, None, None] for r in levels])
    coupling_table = table(["prospective time", "phi L2 rel", "H1 rel", "energy rel", "D_CH rel", "E_kin/excess", "hydroD/D_CH"],
        [[t]+[None]*6 for t in selections["probe_times"]])
    work = audit["rejected_BE_work_diagnostic"]
    work_table = table(["n", "dt", "DeltaF", "Dold", "Dnew", "B_BE", "weak defect", "trapezoid gap", "continuum defect"],
        [["1 REJECTED", work["dt"], work["Delta_E"], work["D_old"], work["D_new"], work["B_BE"],
          work["weak_work_defect"], work["trapezoid_gap"], work["continuum_local_defect"]]])
    norms = audit["algebraic_norms_diagnostic_only"]
    report = f"""# STEP 3A.4 — cost-efficient CH validation checkpoint

**MODEL NOT YET VALIDATED.** Linear architecture succeeds; nonlinear pilot
stops on the FIRST interval. No accepted stiff-bump trajectory, no three-level
qualification and no full-CHNS coupling transfer. All missing metrics below
are explicitly N/A, not linear-reference substitutes.

## 1. Goal

Separate full sparse LINEAR verification, full NONLINEAR isolated CH, and
targeted full NONLINEAR CHNS. Preserve the physical model and scientific gates.
Implementation base: ec01bf3a570cbf2d58c8e16bc133f43f8a7c8baa.
Policies were frozen in STEP3A4_DESIGN.md before new calculations; policy SHA
`{read('policy/policy.json')['policy_sha256']}`. Historical reports/data unchanged.

## 2. Why STEP3A.3 stopped

The reduced-linear passing 1796-step plan implied 4.65 hours level0 and
{old_series_h:.8f} hours for full-CHNS plan/half/quarter, above the unchanged
3/10-hour old limits. That was a cost checkpoint, not nonlinear failure evidence.

## 3. Old cost estimate

The conservative historical maximum was {old_s:.11f} seconds per accepted
full-CHNS interval. It was inappropriately also used to suppress sparse-linear
verification. It is NOT an isolated-CH timing and NOT measured total CPU usage.

## 4. Separation of analysis/isolated/full-CHNS costs

New fixed budgets: analysis 1 h/stage and 2 h total; optimizer 30 min CPU with
finite evaluation/round/split guards; isolated level0 2 h, pilot+three levels
8 h total; targeted CHNS probe 2 h. These roles have separate authorizations.
Actual charged analysis: {analysis_wall:.6f} s. No budget was increased.

## 5. Full-sparse verification independent of nonlinear cost

The old `spectral_ch.planner_stage` put `verify_full_sparse_schedule` inside
`if cost_ok`, where cost_ok meant the FULL CHNS series cost. The new STEP3A4
entry point uses `run_sparse_analysis`: reduced PASS plus its OWN analysis
budget, with no nonlinear cost input. The regression explicitly uses an
unaffordable CHNS series and verifies that the linear callback still executes.
The same dependency is removed from the legacy planner entry point and tested
directly; its legacy authorization refers only to a full-CHNS series. Historical
STEP3A3 artifacts remain unchanged and are not regenerated.
Both mandatory sparse verifications actually executed and PASS all five gates.

## 6. Timing profile

{cost_table}

Baseline: 19 numerical setups for 19 distinct dt; setup/factorization
{timings['baseline_sparse']['categories']['factorization_setup']['wall_s']:.6f} s;
solves {timings['baseline_sparse']['categories']['linear_solve']['wall_s']:.6f} s;
reference/diagnostics {timings['baseline_sparse']['categories']['diagnostics']['wall_s']:.6f} s;
peak RSS {timings['baseline_sparse']['peak_RSS_bytes']/2**20:.3f} MiB.
Optimized: {timings['optimized_sparse']['categories']['factorization_setup']['count']} setups for
{timings['optimized_sparse']['distinct_factors']} distinct dt (adjacent equal-dt blocks
are independently factorized in this verifier). Sparse solves reuse factors
within each block; no dense production L or inverse is formed.

At t=0, mean ordinary diagnostics 0.153 s; exact/certified P2 0.122 s;
geometry including duplicate resolution 0.755 s; BE work 0.137 s. Contact fits
were already designated sparse in DESIGN. Hard mass/energy/power/resolution
remain every accepted state; exact bottom-facet topology adds ~0.011 s at t=0.
No P2 mathematics was simplified. Cached wall geometry/maps reproduce uncached
and full-contour roots in tests. Further vectorization was not necessary here.

PETSc events are inclusive/nested, unlike application timers; they are not
blindly summed. The failed SNESSolve event has zero completed event time after
exception, so the application failure-safe timer ({timings['isolated_failed_pilot']['categories']['nonlinear_SNES']['wall_s']:.6f} s)
is authoritative for that attempt. See the [PETSc profiling manual](https://petsc.org/release/manual/profiling/).

## 7. Old 898/1796 plan

{plan_table}

The 898 schedule is reconstructed by reversing the historical uniform split,
without changing any block boundary. The reduced quadratic D integral uses
all endpoint powers including t=0, never -DeltaE.

## 8. Blockwise optimization algorithm

Each legal single-block split propagates the reduced solution all the way to
1e-4 and evaluates all four global errors. Rank decrease of maximum normalized
constraint violation per added interval; tie-break smaller block index.
Hard acceptance still requires EVERY limit. Adjacent dt growth <=2 and the
original power-of-two hierarchy remain enforced. Deterministic reverse merges
stop at a fixed point. This is NOT a proof of global optimality.

## 9. Optimization history

{history_table}

All nine tested/accepted candidate records are retained. Only blocks 17 and
16 were refined, in that order. The result is a legal one-block merge fixed
point. Independent repeat semantic hashes match. Optimization and repeat took
{timings['optimizer_and_repeat']['wall_s']:.6f} s including setup, well inside guards.
The repeat compares the deterministic optimizer CORE semantic SHA
`{optimization['plan_sha256']}` (fresh search history/cache per invocation).
The frozen execution-artifact SHA below additionally includes archived source
provenance. Gzip archive timestamps are not part of optimizer reproducibility.
All six scientific source archives were independently checked against both
recorded source hashes and archive-file hashes; all match.

## 10. Optimized plan

1068 intervals: 728 fewer, **40.534521% shorter** than 1796. Half/quarter
would be 2136/4272 intervals with identical boundaries. dt0={plan['dt0']:.16g} s.
Plan semantic SHA: `{plan['plan_sha256']}`.

{blocks_table}

## 11. Reduced reference verification

The qualified historical rational M-Krylov reference, matrices and perturbation
are reused, not regenerated or relabelled as newly generated evidence. Inputs
and reference arrays are fingerprint checked. Reduced optimized gates PASS.

## 12. Full-sparse optimized verification

{sparse_table}

All five full-sparse gates PASS. Actual runtime {timings['optimized_sparse']['wall_s']:.6f} s,
not estimated from old CHNS timing. L=-mobility M0^-1 K M0^-1 H is unchanged;
the sparse two-field solve is the exact balanced linear BE system.
Projection remains analysis-only. These are linear results, not nonlinear evidence.
Linear q0=psi is unscaled: energies/powers in raw linear artifacts require a
factor delta^2=1e-6 for physical comparison with the nonlinear bump.

## 13. Isolated CH chemical-equation equivalence

Tiny automated tests independently assemble and compare phase and chemical
blocks of IsolatedCH and CHNSSolver at u=0 with nontrivial phi/mu, theta=60,
P2 and quadrature 12. Relative difference must be <=1e-12. The same ordinary
weak initial-mu projection agrees coefficientwise. Previous accepted phi/mu
are retained as Newton guesses; there is no reset to equilibrium or clipping.
Tiny three-level nonlinear CH and small-amplitude CHNS coupling tests PASS;
they are smoke evidence only, never production qualification.

## 14. Isolated stiff-bump level0

**NOT RUN.** The predeclared 50-interval cost pilot failed at interval 1,
dt={audit['dt']:.16g} s. Accepted intervals=0. SNES reason=-6
(DIVERGED_LINE_SEARCH), 11 iterations; final residual={pilot['SNES_residual']:.12g}.
Initial residual=0.784388568715275. With unchanged atol=1e-10, rtol=1e-9,
stol=0, the residual target is max(atol, rtol*initial)=7.84388568715e-10.
The stored reassembly reproduces the reported residual. No retry, smaller dt,
tolerance change or solution correction was attempted. See the
[SNES tolerance definitions](https://petsc.org/release/manualpages/SNES/SNESSetTolerances/).

At t=0: phase mass={initial['phase_mass']:.16g}, mass/target-domain error
{initial['mass_error_to_target_domain']:.9g}, certified cells
{initial['cells_across_transition_certified_min']:.12g}, two wall crossings.
D_CH(0)={initial['CH_dissipation']:.17g} W/m and excess
{initial['excess_free_energy']:.17g} J/m reproduce historical values exactly.
Initial phi SHA and psi SHA are unchanged. The DOF range is
[{initial['phi_min_dof']:.12g}, {initial['phi_max_dof']:.12g}]; no clipping.

Initial phi SHA: `{input_provenance['initial_phi_sha256']}`.
Psi SHA: `{input_provenance['perturbation_coefficient_sha256']}`.
Prepared equilibrium: `{input_provenance['equilibrium_fingerprint']}`.
Mesh: `{input_provenance['mesh_fingerprint']}`. Quadrature degree 12, P2;
all physical/config parameters are loaded unchanged from the prepared cache.

Pilot wall time {timings['isolated_failed_pilot']['wall_s']:.6f} s, of which failed
SNES attempt {timings['isolated_failed_pilot']['categories']['nonlinear_SNES']['wall_s']:.6f} s.
This is NOT seconds per accepted interval. Isolated throughput and full-series
cost cannot be established from zero accepted steps. Cost authorization remains false.

## 15. Half refinement

NOT RUN, blocked by pilot failure. No alternate startup or dt retry.

## 16. Quarter refinement

NOT RUN. No formal three-level convergence evidence exists.

{series_table}

## 17. Integrated D_CH convergence

Not measured: no accepted isolated intervals. No Richardson extrapolation,
linear D2 integral or rejected-state endpoint power substitutes for the missing
nonlinear physical-dissipation integrals. Last relative D_CH difference: N/A.

## 18. Energy closure

Nonlinear global closure and successive defect decrease: N/A. The initial
energy is {initial['E_total']:.17g} J/m. A near-zero work residual for one
REJECTED Newton iterate (below) does not qualify a nonlinear time history.

## 19. Endpoint phi convergence

N/A: no level0/half/quarter endpoints. No conditional transfer is inferred.

## 20. Early BE work decomposition

There are NO accepted BE-work intervals. The following is a read-only
postmortem of the stored REJECTED iterate, clearly not trajectory evidence:

{work_table}

Decomposition roundoff={work['decomposition_roundoff']:.9g} J/m.
Physical continuum dissipation is D_CH+D_visc+D_slip. B_BE is the signed
NUMERICAL BE work remainder, diagnostic-only, never physical heat.

Residual localization (algebraic norms are solver-conditioning diagnostics,
not scientific stationarity norms): phase={norms['phase_literal_UFL']:.12g};
chemical={norms['chemical_literal_UFL']:.12g}. The time and flux vector norms
are ~0.773630519 and nearly cancel. Evaluating the same time term using
M0*(phi_coeff-old_coeff)/dt differs from the literal UFL evaluation by
{norms['literal_time_minus_coefficient_time']:.12g}. Its complete phase residual
is still {norms['phase_coefficient_difference_sparse']:.12g}; merely reassembling
the rejected coefficients is NOT a demonstrated fix. Stable-time UFL vs sparse
assembly differs by {norms['coefficient_time_UFL_minus_sparse']:.12g}.

The L2 discrepancy between subtracting interpolated fields and interpolating
coefficient differences is {audit['interpolation_subtraction_discrepancy_L2']:.12g};
division by tiny dt amplifies it to {audit['interpolation_subtraction_discrepancy_L2_over_dt']:.12g}.
Phase Riesz L2 residual={audit['phase_Riesz_L2_literal']:.12g}.
The data support a floating-point time-difference/absolute-state representation
conditioning diagnosis; they do not prove that any proposed representation
change would solve it. A separate accuracy-audited algebraic formulation study
(e.g. retained increments) is needed before any new PDE run. No such solver
change or gate relaxation was made in this iteration.

## 21. Isolated-CH temporal verdict

**NOT QUALIFIED — FIRST INTERVAL SNES FAILURE.** This is a new numerical
blocker after linear temporal resolution passed, not a cost-only stop and not
evidence requiring a change of free energy or continuum physics.

## 22. Coupling probe design

Fixed targets: 5*tau_fast, tau_E, 5*tau_E, 1e-5 s. The predeclared snapping
rule gives the prospective times in the table below (indices 100/590/679/745).
They were determined before pilot execution; an actual CHNS execution plan
was NOT authorized/frozen because the isolated prerequisite failed. The
prospective prefix estimate is {probe['historical_cost_prefix_estimate_wall_s']/3600:.6f} h,
below the 2 h probe limit; this does not bypass the scientific prerequisite.
Coupling thresholds remain phi/E 1e-3, D_CH 1e-2 above 1e-8*D_CH(0),
kinetic/excess and hydroD/D_CH 1e-4, mass 1e-10, certified cells 8.

## 23. Full CHNS coupling results

{coupling_table}

NOT RUN, as required after isolated failure. Initial-condition and schedule
equivalence alone do not demonstrate negligible coupling.

## 24. Hydrodynamic power ratios

No new full-CHNS measurements. The isolated observer's u=0, Ekin=0 and
hydrodynamic powers=0 are CONSTRUCTION of the benchmark, not measured CHNS
zeros. Historical small ratios remain historical, not new transfer evidence.

## 25. CHNS-vs-isolated field differences

N/A at every proposed coupling checkpoint. Do not generalize old tiny coupling
to this resolved production path, later contact-line motion or gravity/sloshing.

## 26. Combined uncertainty

Verified linear planning errors are tabulated separately. Isolated nonlinear
temporal error and full-CHNS coupling discrepancy are both unknown. There is
no combined full-CHNS uncertainty bound or conditional transfer qualification.

## 27. Computational cost saved

The optimizer removes 728/1796 intervals per level (40.53%), with measured
full-sparse runtime reduced from {timings['baseline_sparse']['wall_s']:.6f} to
{timings['optimized_sparse']['wall_s']:.6f} s. Hypothetical full-CHNS plan/half/quarter
would drop from {old_series_h:.6f} to {7*1068*old_s/3600:.6f} wall hours, still not authorized.
New measured scientific stages consumed {wall:.6f} wall s and {cpu:.6f} process
CPU s, excluding test/report and container-launch overhead. Stopping avoided
unqualified expensive runs, but **no completed-validation CPU savings** can
be claimed: the old 32.54 h estimate is wall time and the new validation is incomplete.

## 28. Remaining limitations

The nonlinear first-interval solve must be addressed without weakening the
fixed accuracy policy. Three-level isolated convergence and the coupling
bridge remain unperformed; coupling runner/series qualification beyond the
stop point are not claimed complete. No moving-contact, film, tank, gravity,
water-air or unequal-density runs were started. The prepared equilibrium,
physics, mesh, P2 degree, quadrature and initial coefficient hashes are intact.

Local pure suite: {tests['pure']['passed']} passed, {tests['pure']['skipped']} skipped
(FEniCSx dependencies), 2 validation-marked deselected. Pinned Docker workflow:
{tests['pinned_Docker_workflow']['passed']} passed, {tests['pinned_Docker_workflow']['failures']} failures,
{tests['pinned_Docker_workflow']['errors']} errors. The initial Docker invocation
ran no tests due to host glob expansion; corrected invocation is separately
recorded. GitHub workflow includes STEP3A4 tiny tests; remote GitHub CI was
NOT RUN because no push was authorized. JUnit files and per-stage source
archives/hashes are retained. Older artifact provenance is not rewritten to
claim this HEAD generated it.

## 29. Scientific verdict

**MODEL NOT YET VALIDATED.**

LINEAR SCHEDULE VERIFIED; NONLINEAR ISOLATED CH STOPPED AT FIRST SNES
INTERVAL; NO FULL-CHNS TRANSFER QUALIFICATION.

This checkpoint does NOT claim NONLINEAR ISOLATED CH STIFF TRANSIENT QUALIFIED,
BASIC FULL CHNS TRANSIENT QUALIFIED, MODEL VALIDATED or TANK READY.
"""
    Path("STEP3A4_REPORT.md").write_text(report)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for label, blocks in (("uniform 1796", read("linear_sparse/baseline_plan/input_plan.json")["blocks"]), ("optimized 1068", plan["blocks"])):
        x = [max(b["t_start"], 1e-12) for b in blocks]+[blocks[-1]["t_end"]]
        y = [b["dt"] for b in blocks]+[blocks[-1]["dt"]]
        axes[0].step(x, y, where="post", label=label)
    axes[0].set(xscale="log", yscale="log", xlabel="time [s]", ylabel="immutable BE dt [s]")
    axes[0].legend()
    for j, (label, _, errors, _) in enumerate(plans):
        axes[1].bar(np.arange(4)+(j-1)*.25, [errors[k]/v for k, v in LINEAR_LIMITS.items()], width=.25, label=label)
    axes[1].axhline(1., color="black", ls="--")
    axes[1].set(xticks=np.arange(4), xticklabels=["endpoint", "max L2", "energy", "D2"], ylabel="reduced error / HARD limit")
    axes[1].legend(fontsize=8)
    fig.savefig(ROOT/"schedule_comparison.pdf"); plt.close(fig)

    snes = [json.loads(line) for line in (ROOT/"performance/isolated_pilot/snes_iterations.jsonl").read_text().splitlines()]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    axes[0].semilogy([r["iteration"] for r in snes], [r["residual"] for r in snes], "o-")
    axes[0].axhline(7.84388568715e-10, color="red", ls="--", label="unchanged residual target")
    axes[0].set(xlabel="Newton iteration (interval 1)", ylabel="algebraic SNES norm", title="DIVERGED_LINE_SEARCH; zero accepted steps")
    axes[0].legend()
    axes[1].bar(["-DeltaF", "trap integral D", "B_BE", "trap gap"],
        [-work["Delta_E"], work["trapezoid_integral"], work["B_BE"], work["trapezoid_gap"]])
    axes[1].set(yscale="log", ylabel="J/m", title="REJECTED iterate diagnostic ONLY — no trajectory")
    fig.savefig(ROOT/"isolated_energy.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 3), constrained_layout=True)
    ax.scatter(selections["probe_times"], np.ones(4), marker="|", s=300, color="gray")
    for j, t in enumerate(selections["probe_times"]):
        ax.annotate(f"T{j+1}\n{t:.5g} s", (t, 1), xytext=(0, 15), textcoords="offset points", ha="center")
    ax.set(xscale="log", ylim=(.7, 1.4), yticks=[], xlabel="prospective reference-derived checkpoints [s]",
           title="Coupling comparison NOT RUN: isolated prerequisite failed")
    fig.savefig(ROOT/"coupling_comparison.pdf"); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for i, key in enumerate(("baseline_sparse", "optimized_sparse")):
        left = 0.
        for category in ("assembly", "factorization_setup", "linear_solve", "diagnostics", "initialization", "filesystem_write"):
            value = timings[key]["categories"][category]["wall_s"]
            axes[0].barh(i, value, left=left, label=category if i == 0 else None); left += value
    axes[0].set(yticks=[0, 1], yticklabels=["1796 linear", "1068 linear"], xlabel="measured application wall time [s]")
    axes[0].legend(fontsize=7)
    names = list(timings)
    axes[1].barh(names, [timings[n]["wall_s"] for n in names])
    axes[1].set(xlabel="measured wall s (not a completed nonlinear series)")
    fig.savefig(ROOT/"cost_breakdown.pdf"); plt.close(fig)
    print(json.dumps({"verdict": result["scientific_verdict"], "tests": tests,
        "scientific_stage_wall_s": wall, "scientific_stage_CPU_s": cpu}), flush=True)


if __name__ == "__main__":
    main()
