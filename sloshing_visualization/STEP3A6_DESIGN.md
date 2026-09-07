# STEP 3A.6 — fixed one-step accuracy and tiny-step protocol

Base HEAD: 7ed5609cb9f13012f750b3f7fce7c98454386b4c. Overall scientific
verdict remains **MODEL NOT YET VALIDATED**. This design is frozen before
new scientific PDE runs. STEP3A4/5 reports and data are historical/read-only.

## Three distinct accuracies

1. Algebraic equivalence: at corresponding physical trial states, R_rate =
   R_absolute in exact arithmetic; J_rate=J_absolute*diag(dt,I_mu). This is
   a column transformation, not residual-row scaling.
2. Nonlinear root accuracy: an unscaled mixed Euclidean stopping norm does
   not guarantee each field/output error. The primary guesses differ:
   absolute starts at phi_old; rate starts at phi_old+dt*r_semidiscrete.
3. Temporal BE accuracy: neither the moderate comparison nor one tiny step
   establishes convergence of a time trajectory or physical dissipation.

## Immutable physics and input checks

Reuse exact prepared equilibrium, initial weak mu projection, 96x48 mesh,
P2/P2, quadrature 12, delta=1e-3 and the saved psi. No model/free-energy/wall/
material/mobility/slip/forcing change. No clipping, smoothing, mass or energy
correction. No BDF2 or original IsolatedCH/CHNSSolver change. Production
ModelConfig is unchanged; validation SNES controls are separate metadata.

Required hashes, checked against frozen policies AND actual initialized arrays:

- initial phi: 6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681;
- psi: 831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0;
- equilibrium: df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc;
- optimized plan: 5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8.

Reproduce initial D_CH=0.11400313905309394 to relative 1e-10 and certified
cells >=8 (historical 8.328131075834218). Verify two wall crossings.
The configured observer schedule is not the actual one-step experiment.

## Finite list of independent moderate solves

Each run starts fresh at t=0 from exactly the same physical input and mu_old;
each solves ONLY [0,1e-6 s]. No endpoint becomes the next old field.
Historical P0 (atol1e-10/rtol1e-9) is read-only evidence, not rerun.

| level | atol | rtol | absolute guess | rate guess |
|---|---:|---:|---|---|
| P1 | 1e-12 | 1e-12 | phi_old | semidiscrete rate |
| P2 | 3e-13 | 3e-13 | phi_old | semidiscrete rate |
| P2 matched-guess control | 3e-13 | 3e-13 | reuse matching P2 absolute | zero rate |

All controls otherwise fixed: stol=0, max_it=30, newtonls/bt,
preonly/LU/MUMPS. No variable scaling, custom line search or unlisted polish.
Each solve has a unique options prefix and records actual controls before
and after. Effective target is max(actual atol, actual rtol*initial residual).

P1 cross-comparison FAIL permits the predeclared P2 experiments. Nonfinite,
material, mass, resolution, topology or unexplained energy/work failure
stops downstream PDE calculations. A failed P2 Newton attempt remains
REJECTED; accuracy-floor evidence is recorded without accepting its iterate.

Cross-formulation thresholds remain exactly:

- ||phi_rate-phi_abs||_M/||phi_abs||_M <=1e-10;
- ||mu_rate-mu_abs||_M/||mu_abs||_M <=1e-10, without mean subtraction;
- |D_rate/D_abs-1| <=1e-10;
- |F_rate-F_abs| <=1e-12 J/m;
- |mass_rate-mass_abs| <=1e-12 m^2.

For P2 primary and matched-guess solutions, solve J e=-R once diagnostically,
without publishing x+e or advancing time. Report relative physical e_phi
(dt*e_rate in the rate formulation), e_mu, and
deltaD=2*mobility*mu^T K e_mu, quadratic=mobility*e_mu^T K e_mu.
Also report first-variation energy and mass corrections. Engineering targets:
each estimated output error <=one quarter of its comparison threshold
(relative phi/mu/D 2.5e-11; absolute F/mass 2.5e-13). These are local Newton
estimates, not rigorous a posteriori bounds. A sparse solve residual and
finite correction are required; unreliable/floor-limited estimates are
inconclusive, not a guarantee. No ideal order versus tolerance is required;
P1/P2 may perform identical Newton iterations.

Moderate qualification requires algebraic/Jacobian tests, both P2 physical
checks, all cross gates, engineering estimates, no larger unresolved accuracy
change, matched-guess agreement and exact retained-rate archive roundtrip.

## Transactional candidate and archival rules

SNES success -> synchronize actual fem.Function/PETSc iterate -> reconstruct
candidate -> finite/material checks on a separate Function -> physical gates
on an isolated diagnostic observer -> publish accepted state/time/history.
Failures preserve state, old, older, counters, current dt/phase, previous
accepted rate/history and cumulative diagnostics. No retry. Original
absolute weak forms are reused in a one-step validation wrapper, not edited.

Archive independent copies of old phi, actual mixed rate, new mu, rounded new
phi, dt*rate, physical-state difference and recovered-rate DIAGNOSTIC. Include
scalar/mixed DOF maps, mesh/layout/version, SHA256 of every array and actual
PETSc iterate. Successful snapshots require synchronized arrays and residual
reassembly after roundtrip. Rejected Function / SNES iterate / accepted
physical state have distinct names and accepted=false metadata.

## Physical gates (all newly accepted one-step candidates)

Finite arrays/metrics and existing material admissibility; mass/domain <=1e-10;
certified cells >=8; two wall crossings; energy growth <=1e-9 J/m;
|weak-work defect| and |decomposition roundoff| <=1e-12 J/m. Independent
groups: absolute_physical_checks, rate_physical_checks,
cross_formulation_checks, nonlinear_accuracy_checks. A cross PASS cannot
override either physical FAIL. Include t=0 in all physical power accounting.

## Linear tiny comparator and exact target

Only after moderate qualification: one old balanced full-sparse linear BE
and one rate-block solve at dt=4.657657028534679e-12 s. No spectrum or plan
rerun. Historical q0 is mass-projected psi WITHOUT delta; verify that
convention and archive it. Gates: endpoint relative L2<=1e-10, increment
relative L2<=1e-6, mass/domain<=1e-11. Retain dt*rate and both rounded
endpoint differences. No nonlinear execution if the linear comparator fails.

Then ONE primary isolated nonlinear rate step from fresh t=0 at the EXACT
target dt. Restore production controls explicitly: atol=1e-10, rtol=1e-9,
stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS. Primary guess semidiscrete.
Require positive SNES reason, physical gates, exact input hashes and archives.
No changed dt/options or retries on failure.

Only after primary acceptance/physical PASS: one rate-zero control, then
one fresh absolute negative control, same target dt and production controls.
Both start at the historical absolute physical guess. The optional absolute
control is elected in this design (exactly one attempt). Its anticipated
divergence is diagnostic, not a primary failure or a reason to discard
accepted rate evidence. A failed primary blocks all controls. No repeated
negative control to obtain an expected outcome.

## Target read-only proof and interpretation

For each accepted rate archive reassemble, without Newton:
retained-rate UFL/sparse, original literal absolute UFL, coefficient-difference
sparse, chemical on old+dt*retained rate and on rounded physical Function.
Report every block's algebraic and M-Riesz norms. A large literal absolute
residual is diagnostic, not an acceptance veto on the rate formulation.

Record both DeltaF evaluations, physical old/new D, trapezoid/endpoint
integrals, signed B_BE, weak/continuum defects and decomposition. B_BE is
numerical work, never physical heat. Weak-work normalized by interval
activity max(|DeltaF_work|,dt*(Dold+Dnew)/2,1e-30) is diagnostic only.
Report reconstruction increment discrepancy (M-L2/max/relative), ULP/dt,
and Decimal rounding audit of deterministic worst-discrepancy coefficients.

Compare dt*nonlinear_rate to delta*(q_new_linear-q0), delta exactly once.
Keep previous nonlinear relative increment limit 1e-3. Report retained-linear
increment comparison as well as rounded linear subtraction. No stronger
linear/nonlinear equality claim. Compare old rejected absolute only at
matching dt/input hashes, never copy it as an initial guess.

Cancellation support requires primary+zero-rate success, full physical/
linear/reconstruction checks, absolute failure at matching inputs, and
literal phase residual >10*retained phase residual AND actual rate SNES
target. If only semidiscrete guess succeeds, representation-only cause is
not isolated. If absolute also converges, historical failure is not reproduced.

## Budget, sequence and status

Total new scientific wall budget **1800 s (30 min)**, including scientific
initialization/JIT, diagnostic sparse solves and archival I/O; excluding
container build and general pytest. Check before each independent run and
charge all measured stages. No unbounded repeats or extra accuracy levels.
Separate initialization/JIT, assembly, SNES, diagnostics and archival timers;
PETSc inclusive events are reported without double counting.

Sequence: freeze -> transactional/snapshot/control fixes -> cheap tests ->
P1/P2 -> matched guess + correction estimates -> qualified linear tiny ->
primary tiny -> controls/read-only analysis -> report/tests/commit -> STOP.
CLI --all covers only this sequence. No pilot, 1068-step plan, half/quarter,
full CHNS rate solver, coupling probe or moving-contact/film/tank work.

Statuses are separate: execution complete/failed/not_run; solver converged/
diverged/not_run; physical passed/failed/not_run; moderate qualified/failed/
inconclusive/not_run; target accepted/rejected/not_run; hypothesis supported/
inconclusive/not_reproduced. pilot_authorized=false, full_series_authorized=false,
full_CHNS_transfer=not_qualified, always MODEL NOT YET VALIDATED.

Each result records actual current source/HEAD/archive, physical config,
historical provenance in a separate object, immutable fingerprints, actual
controls/guesses, precision/stack, all state hashes and real reasons. Historic
fields must never overwrite current implementation provenance.
