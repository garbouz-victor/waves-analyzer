# STEP 3A.5 — phase-rate conditioning checkpoint

Overall scientific verdict: **MODEL NOT YET VALIDATED**.

**STOP at the production moderate-dt equivalence gate.** Both formulations
converged, but the prescribed mu and D_CH comparison thresholds failed.
The production tiny-dt experiment, pilot and full-CHNS rate experiments were
therefore NOT RUN. This is not evidence that the tiny-dt rate formulation
fails; its production proof was not authorized by the prerequisite gate.

Measured data: [summary.json](validation_results/step3a5/summary.json),
[summary.csv](validation_results/step3a5/summary.csv). Policy was frozen in
[STEP3A5_DESIGN.md](STEP3A5_DESIGN.md) before new PDE runs; its machine record
and immutable design hash are in [policy.json](validation_results/step3a5/policy/policy.json).

## 1. Goal

Test a BE algebraic change of variables without changing physics, FE spaces,
quadrature, solver tolerances, prepared equilibrium or the optimized schedule.
No clipping, smoothing, mass/energy correction, omitted initial power, altered
line search or timestep retries were used. No three-level scientific series
was authorized or executed. The old IsolatedCH, CHNSSolver and BDF2 are unchanged.

## 2. STEP3A4 first-step blocker

Historical rejected absolute step, **not a solution**:

| metric | absolute old | phase-rate new at the same tiny dt |
|---|---:|---|
| dt [s] | 4.657657028534679e-12 | NOT RUN |
| SNES reason | -6, DIVERGED_LINE_SEARCH | NOT RUN |
| Newton iterations | 11 | NOT RUN |
| initial residual | 0.784388568715275 | NOT RUN |
| final residual | 1.110911036727894e-8 | NOT RUN |
| phase algebraic residual | 1.110911036727894e-8 | NOT RUN |
| chemical algebraic residual | 5.841824769844144e-17 | NOT RUN |
| phase Riesz L2 | 2.836425012469171e-6 | NOT RUN |
| SNES wall time [s] | 12.620930377976038 | NOT RUN |
| accepted? | no, zero accepted intervals | no experiment |

Source: [historical failed pilot](validation_results/step3a4/performance/isolated_pilot/status.json)
and [historical rejected-iterate audit](validation_results/step3a4/performance/failed_interval_audit/status.json).
These artifacts and STEP3A4_REPORT.md were not rewritten or relabelled.

## 3. Floating-point cancellation diagnosis

The STEP3A4 hypothesis remains plausible: phi coefficients near unity lose
digits when subtracted and divided by a time interval of order 1e-12 s.
A pure 70-digit Decimal scalar regression demonstrates that loss of accuracy
while the rate residual retains the root. A tiny nonlinear FEM rate smoke
also converges. Neither substitutes for the required 96x48 same-step proof.

The new moderate comparison exposes an earlier accuracy gate: SNES convergence
alone does not guarantee agreement at 1e-10 relative in every physical field.
The read-only postmortem below investigates that mismatch without a new solve.

## 4. Algebraic change of variables

Old unknowns: (phi_new, mu). Old weak residuals:

```
((phi_new - phi_old)/dt, s) + mobility*(grad(mu), grad(s)) = 0
(mu, chi) - F'_h(phi_new)[chi] = 0
```

New unknowns: (phase_rate, mu), both P2. New weak residuals:

```
phi_new_expr = phi_old + dt*phase_rate
(phase_rate, s) + mobility*(grad(mu), grad(s)) = 0
(mu, chi) - F'_h(phi_new_expr)[chi] = 0
```

[phase_rate.py](src/sloshing/multiphase/phase_rate.py) calls the existing
chemical_stationary_form, including the exact current wall derivative and
quadrature 12. The nonlinear residual never subtracts new and old phase.
Rate and mu unknowns are separate from the physical observer. Physical
coefficients are reconstructed only after positive SNES convergence.
Failed SNES (exception or negative reason) publishes nothing and advances
neither time nor accepted count. This is tested, including non-aliasing.

## 5. Mathematical equivalence

For dt>0 the affine map phi_new=phi_old+dt*r is invertible in exact arithmetic.
Substitution into the old equations yields the new equations with no altered
functional or parameter. With matched row/column ordering:

```
R_rate(y) = R_abs(phi_old + dt*r, mu)
J_rate = J_abs * diag(dt*I_phi, I_mu)
```

The phase time Jacobian becomes the mass matrix rather than M/dt. This is
a change of unknowns, not a physical stabilization, modified dissipation or
mass correction. Different floating-point evaluations and Newton starting
states need not terminate at numerically identical approximations.

## 6. FE residual equivalence

Tiny FEM arbitrary-state tests pass assembled residual and Jacobian-chain-rule
comparisons at relative tolerance 1e-12 (dt=0.125), and central-difference
Jacobian validation at 1e-6. These tests use the same chemical form.

On saved production moderate states, chemical residual evaluations from
original and reconstructed rate expressions differ by only 6.32e-17 and
6.27e-17 in algebraic norm. The physical reconstructed coefficients match
bitwise in this diagnostic. No production dense Jacobian was formed.

## 7. Moderate-dt comparison

Production 96x48 mesh, ONE step per formulation, dt=1e-6 s. Same exact
initial coefficients and consistent mu. Primary rate initial guess was the
unprojected semidiscrete CH rate; absolute retained its original initialization.
Both use atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS.

| metric | absolute | rate | difference / comparison |
|---|---:|---:|---:|
| phi L2 | 0.8067374409515264 | 0.8067374409515263 | relative field difference 1.975869479447110e-11; PASS |
| mu L2 | 0.1526713997434452 | 0.1526713992038911 | relative field difference 1.672496974051572e-8; FAIL |
| phase mass | -0.6060961027499434 | -0.6060961027499463 | absolute 2.886579864025407e-15; PASS |
| F [J/m] | 0.06564277766876461 | 0.06564277766876432 | absolute 2.914335439641036e-16; PASS |
| D_CH [W/m] | 0.003669671891309622 | 0.003669671885484852 | relative 1.587272535630291e-9; FAIL |

Field differences above are norms of coefficient-vector differences in the
consistent mass metric, **not differences of the displayed norms**.
Frozen relative thresholds are 1e-10 for phi, mu and D_CH; absolute F and
mass thresholds are 1e-12. No thresholds were changed after measurement.

| Newton metric | absolute | rate |
|---|---:|---:|
| reason | +2, CONVERGED_FNORM_ABS | +2, CONVERGED_FNORM_ABS |
| iterations | 1 | 2 |
| initial residual | 0.784388568715275 | 0.002249909454089338 |
| final residual | 1.768754327738801e-11 | 1.068782824220492e-14 |
| atol | 1e-10 | 1e-10 |
| rtol times initial residual | 7.843885687152750e-10 | 2.249909454089338e-12 |
| effective target | 7.843885687152750e-10 | 1e-10 |
| SNES wall time [s] | 1.432519044989022 | 2.047763762995601 |
| complete step call [s] | 1.464753159001702 | 2.063674739009002 |

Raw bt traces show full steps, lambda=1, without rejected line-search trials.
Each KSP solve reports one preonly iteration, reason=4; MUMPS INFOG(1)=0.
Jacobian infinity norms are 40.18055555555626 (absolute) and
10.66669618055575 (rate), not condition numbers. The combined comparison takes
24.8928 s wall / 15.1324 s CPU including initialization, JIT and diagnostics;
process peak RSS is 757026816 bytes. Timing events are nested and not summed
as independent costs. Detailed [timing](validation_results/step3a5/performance/timing_breakdown.json).

## 8. Linear sparse tiny-step comparison

**NOT RUN**, blocked by section 7. A gated analysis command is implemented
for the original balanced sparse BE and the linear rate block system:

```
[ M          mobility*K ] [v ] = [0   ]
[-dt*H       M          ] [mu]   [H*q0]
q_new = q0 + dt*v
```

Its production errors are not measured and no qualification is claimed.
A regression verifies that a failed moderate status prevents even assembly.

## 9. Production-size tiny-step SNES result

**NOT RUN**. The target remains exactly 4.657657028534679e-12 s; it was not
enlarged, reduced, replaced by the moderate run, or bypassed by a pilot.
The +2 reason in section 7 applies ONLY to dt=1e-6, not the target tiny step.

Measured pre-solve semidiscrete guess, on the identical initial input:
rate L2=283.7841311245012, max|rate|=15786.77034436708, mu L2=0.1526736375081269.
At the moderate dt, dt*max|rate|=0.01578677034436708.
Accepted moderate rate L2=6.065135024042037, max|rate|=75.85024192982911,
dt*max|rate|=7.585024192982910e-5. These are not a tiny-step solution.

## 10. Residual decomposition on same physical state

The required table at a **production tiny accepted state** has no values:

| residual evaluator | tiny accepted-state norm |
|---|---|
| retained rate formulation | NOT RUN |
| original literal UFL subtraction | NOT RUN |
| coefficient-difference / dt | NOT RUN |
| sparse mass-matrix | NOT RUN |
| phase/chemical Riesz L2 | NOT RUN |

Instead, a labelled **read-only moderate postmortem** reassembles both saved
physical solutions without any Newton solve, time advance or correction:

| evaluator / metric | saved absolute state | saved rate physical state |
|---|---:|---:|
| original literal phase algebraic | 1.139524294204147e-13 | 1.137822323619331e-13 |
| coefficient difference UFL phase | 8.810893774038857e-14 | 8.834685804640059e-14 |
| recovered rate expression phase | 8.810907947909439e-14 | 8.834694721991113e-14 |
| sparse coefficient phase | 1.012258393667909e-13 | 1.017906875256587e-13 |
| physical chemical algebraic | 1.768717620292462e-11 | 3.551209351639062e-15 |
| literal phase Riesz L2 | 2.955479797368153e-11 | 2.966794974788412e-11 |
| coefficient phase Riesz L2 | 2.204999561939573e-11 | 2.218840452846809e-11 |
| chemical Riesz L2 | 3.644600265229293e-9 | 7.443809951244999e-13 |

Riesz norms are sqrt(R^T M^-1 R), not Euclidean SNES norms. Moderate rate
coefficients were not archived: the postmortem rate is explicitly recovered
as (saved_phi-initial_phi)/dt for **diagnosis only**, never assigned to a
physical solution. Therefore its residual is not the retained unknown's
SNES residual; the latter's measured combined norm remains 1.068782824220492e-14.

One sparse Jacobian **matvec**, not a Newton solve, checks the measured
displacement. Chemical R_abs+J_abs*(x_rate-x_abs)-R_rate has relative norm
2.474431876722029e-6 compared with R_abs, i.e. near assembly roundoff in
absolute units. The mu difference L2 is 2.553424540951301e-9; its mean is
-6.358996293444232e-10 and zero-mean L2 is 2.495762745643122e-9.

Inference: the moderate discrepancy is consistent with different attained
Newton stopping accuracies. Absolute stops after one iteration with a small
but nonzero chemical defect; rate requires two iterations from its different
initial guess. This supports an accuracy-audit explanation, **not** permission
to overrule the fixed comparison gate or a proof of the tiny-dt root cause.
Full [read-only audit](validation_results/step3a5/algebraic_equivalence/moderate_postmortem/status.json).

## 11. ULP amplification analysis

At the saved moderate states max coefficient ULP/dt=2.220446049250313e-10.
Historical tiny rejected state: 4.767302606540087e-5. These coefficient-scale
figures are not directly comparable to assembled Euclidean residual norms.
The historical interpolation/subtraction discrepancy divided by dt was
7.829962076753288e-6 in L2; its phase Riesz defect was 2.836425012469171e-6.
No new accepted-tiny-state ULP measurement exists.

## 12. Energy/work result

Only the moderate accepted states are newly measured here. Both have two
wall crossings and 8.328131075834218 certified transition cells. Mass change
over domain is 7.40149e-15 (absolute) / 3.39235e-15 (rate). The initial D_CH
reproduces 0.11400313905309394 exactly. Initial and final physical fields are
finite, without clipping. Free energy decreases in both cases.

| moderate one-step metric | absolute | rate |
|---|---:|---:|
| DeltaF, subtract saved total energies [J/m] | -5.015781665007602e-9 | -5.015781956441145e-9 |
| DeltaF, assembled BE-work difference [J/m] | -5.015781591774159e-9 | -5.015781588651945e-9 |
| D_old [W/m] | 0.11400313905309394 | 0.11400313905309394 |
| D_new [W/m] | 0.003669671891309622 | 0.003669671885484852 |
| trapezoid integral [J/m] | 5.883640547220178e-8 | 5.883640546928939e-8 |
| endpoint BE integral [J/m] | 3.669671891309621e-9 | 3.669671885484852e-9 |
| B_BE [J/m] | 1.346109707294476e-9 | 1.346109703153599e-9 |
| weak-work defect [J/m] | 6.811477897481704e-18 | -1.161568375177589e-20 |
| trapezoid gap [J/m] | 5.516673358089216e-8 | 5.516673358380454e-8 |
| continuum local defect, BE-work evaluator [J/m] | 5.382062388042762e-8 | 5.382062388063745e-8 |
| decomposition roundoff [J/m] | 1.845605382728315e-20 | -1.872736906820055e-21 |

The two DeltaF evaluations are reported separately, not silently substituted.
Physical continuum dissipation remains D_CH+D_visc+D_slip. Here velocity,
kinetic energy, D_visc and D_slip are zero **by isolated-CH construction**,
not measured full-CHNS negligible coupling. B_BE is a signed numerical BE
work remainder, diagnostic only, never physical heat. The large moderate
trapezoid defect is not a temporal qualification; this is one equivalence step.

## 13. Comparison with historical rejected iterate

No same-dt accepted rate solution exists, so field differences to that iterate
are **NOT MEASURED**. The following is historical reference only:

| physics metric | old REJECTED tiny iterate | new accepted tiny rate state |
|---|---:|---|
| DeltaF from BE-work evaluator [J/m] | -5.235196457849968e-13 | NOT RUN |
| D_new [W/m] | 0.11187268676030823 | NOT RUN |
| phase mass | -0.6060961027499477 | NOT RUN |
| phi min/max DOFs | -0.9790806988061016 / 1.0156953397485375 | NOT RUN |
| mu std [Pa] | 0.001642366867567342 | NOT RUN |
| B_BE [J/m] | 2.455017106619540e-15 | NOT RUN |
| weak-work defect [J/m] | -2.794557704908750e-20 | NOT RUN |
| continuum defect, BE-work evaluator [J/m] | 2.506418052936967e-15 | NOT RUN |

These old coefficients were neither copied into a new solution nor used as
a convergence target. A moderate-dt solution is not a valid same-time comparator.

## 14. Pilot result if run

**NOT RUN**, unauthorized after moderate FAIL. No 50-step accepted cost,
pilot cost gate, block-change production evidence or temporal convergence
claim exists. The implemented next-rate guess policy is previous accepted
rate, unscaled at dt changes; tested only by tiny equilibrium/kinematic tests.
The unchanged 1068-step schedule was not executed or edited. No level0,
half or quarter series was launched, regardless of apparent moderate cost.

## 15. Full-CHNS rate formulation

**NOT IMPLEMENTED in this checkpoint**: phase 10 is downstream of isolated
proof/pilot. The design records the scoped matched-density BE alternate
(u,pi,phase_rate,mu), but no untested new full-CHNS implementation is claimed.
Existing production CHNS and BDF2 remain unchanged. A separate continuation
must first handle the moderate accuracy gate under a newly reviewed design.

## 16. Full-CHNS one-step result

Moderate full-CHNS comparison and production tiny full-CHNS rate smoke:
**NOT RUN**. Existing tiny full-CHNS regression tests remain green, but they
do not test a new rate formulation and do not qualify the coupling bridge.

## 17. Remaining limitations

No production same-dt proof, retained tiny-rate residual reassembly, tiny
linear update comparison, rejected-iterate match, pilot or CHNS transfer was
obtained. ABSOLUTE-STATE CANCELLATION ROOT CAUSE SUPPORTED is **not** claimed.
The moderate stopping-accuracy mismatch needs a separately specified accuracy
audit; no extra Newton step, changed guess, tolerances or relaxed criterion
was tried after the failure. No evidence is being extrapolated to tank,
contact motion, gravity, film or unequal density.

Tests at the final code snapshot:

- Outside Docker: **270 passed, 21 skipped, 2 deselected**. The skips require
  FEniCSx/SLEPc; two existing validation-marked cases retain normal CI exclusion.
- Pinned Docker STEP3–STEP3A5 tiny suite: **176 passed**. New STEP3A5: 10 tests,
  comprising 7 pure and 3 FEM checks; no production mesh in CI.
- Initial full suites each exposed one test-only regex mismatch (`failure`
  versus `failed`). The assertion was fixed; original XML and final XML are
  both retained. No solver or scientific threshold changed for that repair.
- [Pure final XML](validation_results/step3a5/pure_tests_final.xml),
  [Docker final XML](validation_results/step3a5/docker_tests_final.xml).
- GitHub workflow now includes test_step3a5_*.py. **Remote GitHub Actions NOT
  RUN**: no push was performed. Local suite results are not remote-CI claims.

All new scientific measurements were generated at actual HEAD
6106f8977be23a12446c2ea2d5aed0ea3311cf17 with separately archived dirty-source
snapshots. The moderate source hash is
870d6df3b8f50b107782732ea9af1f4dca4da5a66fa88f84d3a42ea4c2b6ec35;
postmortem and report have their own provenance. The later checkpoint commit
does not retroactively relabel the measurement source. Historical artifacts
retain their own original provenance.

The `dt_schedule` inherited in physical-observer provenance describes its
configured historical schedule, not execution in this one-off experiment.
Actual moderate dt=1e-6 is recorded in each step log and in the frozen policy;
the optimized plan was hash-checked, **not executed**.

Pinned Docker image:
sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8;
dolfinx 0.10.0, PETSc 3.24.0, MUMPS 5.8.1, UFL 2025.2.0, quadrature 12.
Initial phi, psi, prepared equilibrium and optimized plan SHA256 all match
the four mandatory hashes in the unchanged DESIGN/policy record. Prepared
cache algorithm files, original IsolatedCH and historical results are unchanged.

## 18. Scientific verdict

**MODERATE-dt EQUIVALENCE GATE FAILED; PRODUCTION TINY-dt PHASE-RATE PROOF NOT RUN.**

**MODEL NOT YET VALIDATED.**

No claim that the tiny-dt blocker is resolved, that rate formulation failed
at the target dt, or that a nonlinear transient/pilot/full CHNS is qualified.
This checkpoint preserves the implementation, measured prerequisite failure
and read-only diagnosis for the next explicitly reviewed accuracy experiment.
