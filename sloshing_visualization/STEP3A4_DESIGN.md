# STEP 3A.4 — policy fixed before new PDE runs

Implementation base: `ec01bf3a570cbf2d58c8e16bc133f43f8a7c8baa`.
Overall verdict remains **MODEL NOT YET VALIDATED**. STEP3/3A1/3A2/3A3 data
and reports are historical/read-only. This document and its machine-readable
policy are saved before scientific runs. No later threshold relaxation.

## Invariants

Reuse the exact STEP3A2 prepared equilibrium, STEP3A3 frozen psi/initial phi,
96x48 crossed mesh, P2 phi/mu, degree-12 quadrature, delta=1e-3, unchanged
free energy, wetting wall, mobility and full CHNS/AGG/pressure/slip conventions.
Require old and new input hashes to match coefficientwise, not merely in norm.
No clipping, smoothing, mass/energy correction, time-zero or interval omission.
Isolated CH is a labelled nonlinear benchmark, never a production model change.

## Three separate authorizations and budgets

* Linear analysis: <=3600 wall seconds per bounded stage, <=7200 wall seconds
  total new scientific linear analysis. Optimizer additionally <=1800 CPU
  seconds, <=5000 candidate evaluations, <=24 greedy rounds, <=3 extra splits
  per original block. The fixed finite guards never relax accuracy constraints.
* Nonlinear isolated CH: level0 <=7200 wall seconds; level0/half/quarter plus
  a separate cost pilot <=28800 wall seconds total. Pilot is the first 50
  exact level0 intervals, with no schedule changes or retries. Extrapolation
  uses the larger of pilot mean and historical STEP3A3 isolated-low-mode mean
  SNES time, plus measured mandatory diagnostics and measured initialization.
  The pilot is NOT scientific series evidence. All three full histories start
  independently from the same frozen initial field if cost permits.
* Targeted full CHNS probe: <=7200 wall seconds total. Predict before execution
  from exact prefix interval count and conservative historical 9.3186596433
  seconds/accepted full-CHNS interval; replace by a measured profile only when
  such an independent profile exists, never by an assumed isolated-CH speedup.

Do not raise the prior 3/10-hour full-series budget. Instead, no full-CHNS
three-level series is authorized here. Report wall time AND process CPU time;
historical 32.5428302877 hours is a wall-cost estimate, not measured CPU usage.
Do not equate those measures when reporting savings.

Statuses are distinct: `reduced_linear_verified`,
`full_sparse_linear_verified`, `isolated_CH_execution_authorized`,
`full_CHNS_probe_authorized`, `full_model_temporal_qualification`.
Full sparse verification has no dependency on nonlinear timing or budget.

## Linear verification and profiling

First verify the untouched STEP3A3 1796-interval candidate on the full sparse
operator using the same coupled balanced MUMPS solve and rational reference.
No dense full L or inverse. Mass projection is analysis-only, as in STEP3A3.
Require at every accepted linear checkpoint:

* endpoint relative L2 <=1e-3;
* maximum checkpoint relative L2 <=1e-3;
* maximum quadratic-energy relative error <=1e-3;
* integrated D2 relative error <=1e-2 (trapezoids INCLUDING t=0);
* max mass/domain <=1e-11.

Do not infer integrated D2 from -DeltaE. Compare reduced and full sparse
results and preserve the complete history. Failure stops downstream stages.
Measure sparse block assembly/setup, numerical factorization, solve time,
reference/diagnostic time, output time, peak RSS and distinct factors.

Before optimization, profile initial-state diagnostics without time advance.
Record mass/energy/power, certified resolution, geometry, BE-work and writing
separately. PETSc event times are inclusive/nested and are NOT blindly added
to exclusive application wall timers. Record unavailable categories as null,
not invented zero costs. Source reference: PETSc profiling manual
https://petsc.org/release/manual/profiling/ ; actual pinned API is tested locally.

## Deterministic blockwise optimization

Only after baseline full sparse PASS, reconstruct the original 898-interval
candidate by reversing the archived uniform halving (same block boundaries).
For each current candidate evaluate all four GLOBAL errors after propagation
to t_end=1e-4. Keep ALL hard constraints, not a weighted soft substitute.

For each legal single-block dt/2 split, repeat full reduced propagation and
evaluate violation `max(0,max(error_i/limit_i)-1)`. Rank reductions in this
violation per added interval; exact ties go to smaller block index. If no
move reduces violation, stop this bounded search and retain the already
full-sparse-verified 1796 baseline. Allow only power-of-two descendants of each
original dt, with unchanged boundaries. A split is legal only if adjacent
step growth remains <=2; reductions in dt at transitions are allowed.

After first PASS, attempt reversions of refined blocks in increasing index
order; accept only if every GLOBAL gate stays PASS and growth <=2. Repeat to
fixed point. This is local minimality under these legal moves, NOT a global
optimum. Preserve every tested candidate with split/merge, errors, pass/fail,
step count and deterministic semantic hashes. Timestamps/timing stay outside
the semantic plan. A finite guard retains a best-found passing plan if any.

Adopt the optimized plan only if it uses fewer than 1796 intervals and passes
both reduced and full sparse verification. Otherwise use the old full-sparse
PASS baseline. A reduced PASS/full-sparse FAIL is STOP, never a tolerance edit.
Run optimizer twice on identical inputs and require identical semantic SHA.

## Nonlinear isolated stiff-bump policy

Use full quartic bulk and current cubic wall energy, not the linear operator.
Automated tiny tests compare phase/chemical residual blocks with full CHNS at
u=0, matching spaces, signs, degree-12 quadrature and ordinary weak initial mu.
Target relative residual-block difference <=1e-12 at a representable test dt.
Previous accepted phi/mu are the next Newton guess. SNES tolerances unchanged:
config atol=1e-10, rtol=1e-9, max_it=30, stol=0. No smaller-dt retry. A failure
at tiny dt triggers a residual/conditioning audit, not a silent tolerance,
representation, timestep or physical-model change during a run.

Levels: frozen plan, dt/2, dt/4, identical block boundaries and t_end=1e-4 s.
Formal qualification requires all three levels. Two-level evidence alone is
NOT qualified. No Richardson extrapolation replaces a measured third level.

Every accepted state (including t=0): finite fields/data, phase mass, physical
free energy and D_CH, DOF extrema, exact P2/certified normal-width resolution.
Require mass/domain <=1e-10 and >=8 certified transition cells at EVERY state.
Physical free-energy growth may not exceed the existing 1e-9 J/m allowance.
u=0, E_kin=0 and hydro powers=0 are CONSTRUCTION of isolated CH, not measured
evidence about CHNS. Final closure |DeltaF+integral D_CH|/|DeltaF| <=5%.

Series: absolute final budget defects decrease, successive endpoint L2
differences decrease, successive D_CH integral differences decrease, and the
last relative D_CH integral difference <=5% (normalize by finest integral).
No threshold change after results. Any isolated scientific failure stops the
full CHNS probe. Linear/low-mode evidence cannot substitute for this series.

## Diagnostics policy fixed before profiling/runs

Mass, energy, power, finite/material checks and certified P2 resolution remain
EVERY accepted state. Contact-wall P2 zero-crossing count is also checked at
every accepted state; require exactly two crossings for this sessile drop.
This explicit topology check prevents assuming preservation from a bulk bump.
Apparent contact-angle fits are diagnostic only here: measure at t=0, block
ends, final time and prescribed coupling times, not every small startup step.
This frequency is selected BEFORE profiling; it is never changed mid-run.

BE work: intervals 1..20 plus unique rounded logarithmic indices from 21 to N
(40 log-spaced indices), every block end and the final interval. Preserve
DeltaF, Dold/new, trapezoid and BE endpoint integrals, signed B_BE, weak-work
defect and decomposition roundoff. Require |weak work| <=1e-12 J/m and
|decomposition roundoff| <=1e-12 J/m. B_BE is numerical work, NOT physical heat.
All continuum power integration still uses every accepted time, including t=0.

Write compact scalar histories every state, full worst-cell certification
witnesses every state, and deterministic field checkpoints for endpoints and
probe times. Cache only mesh-invariant transforms/geometry/maps if profiling
warrants it; require cached/uncached equality to roundoff in automated tests.
No mathematical P2 range or width criterion is weakened.

## Conditional full-CHNS coupling bridge

Only after three-level isolated PASS, freeze a probe JSON before CHNS stepping.
Run ONE continuous full-CHNS trajectory from the identical t=0 state with the
optimized LEVEL0 prefix, not independent restarts at late probe times.
Snapshot labels early/mid/late refer to this same trajectory.

Deterministic target times: T1=5*tau_fast_Ritz, T2=tau_E_initial,
T3=5*tau_E_initial, T4=1e-5 s. For T1..T3 choose the first level0 accepted
time >= target; for T4 choose the last accepted time <=1e-5. Require distinct
increasing times. These rules use saved linear/reference inputs only; the
actual immutable times are written before any full-CHNS probe. Save isolated
level0 fields at these indices in advance for exact-time comparison.

At each common time, measure full-minus-isolated L2/H1 differences, free-energy
excess, D_CH, mass, contact geometry, full speed, kinetic energy, viscous/slip
powers. Require relative phi L2 <=1e-3 with denominator
max(||phi_CH-phi_eq||_L2,1e-12*||delta psi||_L2); H1 analog is diagnostic.
Relative excess-free-energy difference <=1e-3 with denominator
max(|F_CH-F_eq|,1e-12*initial_excess). Physical CH power difference <=1e-2
when isolated D_CH>=1e-8*D_CH(0). Below this floor, report absolute differences.
Require E_kin/initial_excess<=1e-4 and (D_visc+D_slip)/D_CH<=1e-4 wherever
full D_CH>=1e-8*D_CH(0). Check these hydro ratios, mass and certified resolution
at EVERY full-CHNS accepted state, not merely selected comparisons. Same mass
1e-10, resolution 8, finite and 1e-9 energy-growth gates. Do not tighten these
limits after seeing tiny errors. If probe fails, no transfer qualification.

Separate temporal isolated-CH error from measured coupling discrepancy.
Any conditional transfer is limited to this small bulk stiff bump and the
measured EARLY coupling window; it is not full-CHNS temporal refinement over
the entire 1e-4 horizon or evidence for general contact-line/sloshing coupling.

## Sequence, provenance and verdicts

Design -> independent cost gates/profile -> baseline full sparse -> optimizer
and reproducibility -> optimized full sparse -> isolated pilot -> isolated
level0/half/quarter -> isolated series verdict -> frozen probe -> one targeted
full-CHNS prefix -> comparison -> report. STOP on any required gate failure.
No moving-contact spreading, film, tank, gravity or unequal density experiments.

Every NEW result uses actual current HEAD and source hash/archive, while copied
historical provenance remains labelled historical. Store equilibrium/input/
operator/Krylov/mesh/quadrature hashes; plans additionally include optimizer,
policy, history and semantic SHA. Separate execution authorizations by role.

At most: NONLINEAR ISOLATED CH STIFF TRANSIENT QUALIFIED; FULL CHNS TRANSFER
NOT YET QUALIFIED if only isolated PASS. If both PASS: BASIC NONZERO
CH-DOMINATED FULL-CHNS TRANSIENT CONDITIONALLY QUALIFIED THROUGH ISOLATED-CH
TEMPORAL CONVERGENCE AND DIRECT COUPLING-NEGLIGIBILITY CHECK; MOVING-CONTACT
STEP 3A BENCHMARKS STILL REQUIRED. Overall model is NOT YET VALIDATED.
