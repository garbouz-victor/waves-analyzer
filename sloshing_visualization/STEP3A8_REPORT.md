# STEP3A8 — Full-CHNS phase-rate bridge checkpoint

## 1. Goal

Measure full-CHNS hydrodynamic coupling independently of the qualified isolated-CH temporal error, without changing the physical model or BE clock.

**Outcome: mandatory STOP at the first production moderate candidate.** The absolute P1 Newton solve converged, but the candidate exceeded the strong-divergence non-regression envelope fixed in STEP3A8_DESIGN before production. It was not accepted. No production rate solve, target tiny step, pilot or full L0 trajectory was run.

## 2. Inherited STEP3A7 qualification

The three isolated histories remain read-only: L0/L1/L2 = 1068/2136/4272 COMPLETE intervals, 7476 accepted primary steps, no rejected intervals. COMPLETE identities, scalar hash chains, schedules and common field/final-checkpoint hashes were reverified. No isolated PDE, optimizer or spectrum was rerun.

| Isolated L0-to-L2 uncertainty | Measured historical difference |
|---|---:|
| Maximum common-time phi L2 | 3.691188008762449e-8 |
| Final integrated D_CH | 7.033312454099817e-11 J/m |
| Final CH free energy | 6.0652038946784614e-12 J/m |

These are potential transfer limits, not measured full-CHNS coupling.

## 3. Remaining isolated-to-full gap

The required continuous full L0-minus-isolated L0 comparison is NOT RUN. No transfer follows from the small velocity of a rejected moderate candidate. Overall MODEL NOT YET VALIDATED.

## 4. Full CHNS phase-rate derivation

Physical layout remains (u, pi, phi, mu): vector P2 / P1 / P2 / P2. A separate private validation unknown is (u, pi, r_phi, mu), with phi = phi_old + dt*r_phi.

For matched density, drho=0 and J=0. Both formulations use the same momentum weak equation:

```text
(rho*(u-u_old)/dt, v)
+ 1/2 [(grad(u)*(rho*u),v) - (grad(v)*(rho*u),u)]
+ (2*eta(phi)*symgrad(u),symgrad(v))
- (pi,div(v)) + (phi*grad(mu),v)
+ sum_Navier (eta(phi)/slip_length*u_t,v_t)_wall = 0.
```

Continuity is (q,div(u))=0. The absolute CH phase row is

```text
((phi-phi_old)/dt,s) - (phi*u,grad(s))
+ mobility*(grad(mu),grad(s)) = 0.
```

The historical UFL implementation evaluates its BE temporal term as a0*phi+a1*phi_old, with a0=1/dt, a1=-1/dt. The new phase row instead evaluates

```text
(r_phi,s) - ((phi_old+dt*r_phi)*u,grad(s))
+ mobility*(grad(mu),grad(s)) = 0.
```

Chemical equation in both: (mu,chi)-F'_h(phi)[chi]=0, including identical bulk, gradient and wetting-wall variation. The generic drho*potential term remains in the copied equation and vanishes only under the explicit config restriction.

No original CHNSSolver, weak_form.residual, free energy, material helper, BDF2 or historical solver was changed. The validation path rejects unmatched density, gravity, acceleration, non-BE and wrong P2/quad12 configurations.

## 5. Algebraic equivalence

Tiny 4x4 arbitrary, nontrivial fields include u!=0, r!=0 and nonconstant mu, theta=60. These are preflight fields, not the production perturbation.

| Residual term | Relative assembled difference | Frozen limit |
|---|---:|---:|
| Momentum | 4.5891081240737357e-17 | 1e-11 |
| Continuity | 0 | 1e-11 |
| Phase | 1.4284101706339971e-16 | 1e-11 |
| Chemical including wall | 1.2777957859433729e-15 | 1e-11 |
| Chemical wall separately | 1.8841558620952612e-17 | 1e-11 |
| Navier separately | 0 | 1e-11 |

PASS on the tiny mesh. This alone does not qualify production root accuracy or transfer.

## 6. Jacobian transformation

At corresponding physical states, residual rows agree in exact arithmetic. The variable transformation changes columns, not row scaling:

```text
J_rate = J_absolute * diag(I_u, I_pi, dt*I_phi, I_mu).
```

Jacobian action relative difference: 4.2729814365672955e-17 (limit 1e-11).

| Central FD h | Relative action error |
|---|---:|
| 1e-2 | 3.9321600354716027e-13 |
| 5e-3 | 9.830769064527317e-14 |
| 2.5e-3 | 2.461479086787062e-14 |
| 1.25e-3 | 6.645935551186878e-15 |
| 1e-5 | 2.83560417973595e-13 |

Second-order truncation decreases to roundoff; final error <1e-7. PASS.

## 7. BC equivalence

The exact original BC objects are passed to the new problem: identical essential velocity constraints, exactly one pressure gauge and no rate/chemical constraint. DOLFINx 0.10 NonlinearProblem stores BCs in its assembly callbacks, not a public .bcs property; preflight checked the explicitly retained supplied BC list and constrained DOF maps.

Material and physical diagnostics run on a private reconstructed candidate before accepted publication. The true mixed rate is copied directly; it is never replaced by a recovered (rounded_phi_new-phi_old)/dt. Tiny archive/restart tests verify this.

## 8. Moderate full-CHNS comparison

Frozen finite series: independent [0,1e-6] solves, P1 atol=rtol=1e-12 and P2 atol=rtol=3e-13, both mandatory unless a physical gate stops progression. No solution sharing. Phi/mu/D_CH comparison limit 1e-10, energy/mass absolute 1e-12. U_ref=sqrt(sigma/(rho*epsilon))=2 m/s and P_ref=sigma/epsilon=4 Pa; scaled velocity/pressure comparison limits 1e-9.

Only **P1 absolute** ran. Actual before/after getters: atol=rtol=1e-12, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS. No controls changed.

| Metric | Absolute P1 (REJECTED physical candidate) | Phase-rate | Difference |
|---|---:|---|---|
| dt | 1e-6 s | NOT RUN | NOT RUN |
| SNES reason | 2, CONVERGED_FNORM_ABS | NOT RUN | NOT RUN |
| Newton iterations | 2 | NOT RUN | NOT RUN |
| Initial norm | 0.7843890864759492 | NOT RUN | NOT RUN |
| Final norm | 1.0817199474444568e-13 | NOT RUN | NOT RUN |
| Effective target | 1e-12 | NOT RUN | NOT RUN |
| u L2 | 7.563407043299337e-9 | NOT RUN | NOT RUN |
| pi L2 | 0.001091516415835943 | NOT RUN | NOT RUN |
| phi L2 | 0.8067374409515263 | NOT RUN | NOT RUN |
| mu L2 | 0.15267139920381542 | NOT RUN | NOT RUN |
| Phase mass | -0.6060961027499453 m² | NOT RUN | NOT RUN |
| Total energy | 0.06564277766876435 J/m | NOT RUN | NOT RUN |
| D_CH | 0.0036696717867036074 W/m | NOT RUN | NOT RUN |
| D_visc | 3.100629956378667e-11 W/m | NOT RUN | NOT RUN |
| D_slip | 2.440535248219972e-23 W/m | NOT RUN | NOT RUN |
| Accepted? | NO | NOT RUN | NOT RUN |

Newton norms: 0.7843890864759492 → 4.056172190337836e-8 → 1.0817199474444568e-13. Each Newton KSP: one iteration, reason 4. This is successful nonlinear convergence, not an accepted scientific step.

P1 rate, P2 absolute and P2 rate are NOT RUN because the predeclared physical/non-regression gate failed before cross-comparison. The production formulation comparison is NOT QUALIFIED, not evidence of mathematical non-equivalence.

## 9. Exact full target tiny step

Full rate dt=4.657657028534679e-12: **NOT RUN**. P0 controls remain configured but were not applied in any production full-rate solve: atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS.

| Quantity | Full CHNS rate | Isolated CH rate |
|---|---|---|
| Tiny convergence | NOT RUN | Historical STEP3A6 PASS |
| phi/mu differences | NOT RUN | Reference only |
| Max speed / E_kin | NOT RUN | Zero by isolated model |
| D_CH / D_visc / D_slip | NOT RUN | Historical, no new comparison |
| Mass / resolution / work | NOT RUN | Historical STEP3A6 PASS |

The optional absolute full tiny control was not elected in the frozen design and was not run. Its failure is not a transfer prerequisite.

## 10. Pilot

NOT RUN; zero accepted production trajectory intervals. No accepted-step cost, stability or drift claim.

## 11. First dt boundary

NOT RUN. The verified schedule has its first real block end at parent interval 122; the planned pilot endpoint is derived as max(50,122+5)=127. These are schedule facts, not observed coupled behavior.

## 12. Restart

Production fork 25→28: NOT RUN. Tiny coupled restart in a new process across a dt change: PASS. Tiny SNES, material-after-convergence and physical-after-convergence failures preserve physical state, old/older, time, accepted count, retained rate, cumulative diagnostics and checkpoint.

Read-only post-STOP reassembly of the saved moderate candidate confirmed strong divergence and weak continuity exactly. The saved Function equals the synchronized PETSc iterate. The accepted state is distinct from this candidate and retains the exact initial phi hash, u=pi=0. No initialized/advanced/Newton solve was used for this reassembly.

## 13. Cost gate

Frozen before production: preliminary moderate/tiny/pilot/fork <=1800 s, continuous full L0 including reused prefix <=14400 s, total <=16200 s. This explicitly predeclared continuous-horizon budget replaces the old short-probe 7200 s budget; no runtime increase occurred.

Scientific wall spent: 19.65485542299575 s. Cost was **not** the stop reason. Pilot forecast, full L0 forecast and actual full L0 runtime: NOT RUN. Failed/rejected candidate timing is not an accepted trajectory seconds-per-step measurement.

## 14. Continuous full L0 trajectory

NOT RUN. Schedule count 1068, final endpoint 1e-4, SHA
`21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31`.
No COMPLETE marker exists; no partial path is counted as complete.

| Metric | Full CHNS L0 | Isolated CH L0 | Isolated CH L2 |
|---|---|---:|---:|
| Final CH free energy | NOT RUN | 0.06564273942723216 | 0.06564273942116695 |
| Integrated D_CH | NOT RUN | 4.3358983982305045e-8 | 4.3288650857764046e-8 |
| Final phi comparison | NOT RUN | Reference | Reference |
| Full-minus-isolated mass | NOT RUN | Reference | Reference |
| Final budget defect | NOT RUN | 1.0166986272259831e-10 | 2.5271534286921684e-11 |

Kinetic/hydro entries for a full trajectory are NOT RUN, not zero.

## 15. Full physical/mass/resolution diagnostics

The following are exclusively the **REJECTED MODERATE CANDIDATE**, not trajectory maxima:

| Metric | Measured | Frozen criterion | Result |
|---|---:|---:|---|
| Strong divergence L2 | 0.000003814612027895975 | <=3.6880704671415807e-6 | FAIL, +3.431104743843716% |
| Relative strong divergence | 0.9403747803293864 | Diagnostic; historical max 0.9301196221025294 | Diagnostic |
| Weak continuity algebraic norm | 4.594698458802348e-24 | Diagnostic | — |
| Weak continuity Riesz L2 | 1.1783432940661996e-21 | <=1e-12 | PASS |
| Mass change/domain | 4.780126911580535e-15 | <=1e-10 | PASS |
| Certified transition cells | 8.328131075834218 | >=8 | PASS |
| Wall crossings | 2 | exactly 2 | PASS |
| Essential velocity / gauge | 0 / 0 | scaled <=1e-12 | PASS |
| Finite / material | PASS | required | PASS |

The failed limit was a STEP3A8 **empirical no-regression envelope** taken from committed same-mesh STEP3A2 stiff-bump histories, not a previously universal strong-divergence tolerance. Its maximum comes from [the dt=2.5e-6 history](validation_results/step3a2/equilibrium_perturbation/be_dt_2_5e-6/history.csv). Those histories have different dt/time samples. Its failure does not contradict the extremely small weak continuity residual and does not by itself establish a continuum-physics error. The frozen protocol nevertheless requires STOP; the envelope was not relaxed or relabelled after measurement. Any reassessment belongs to a separately authorized audit.

## 16. Full energy/work

No continuous full-trajectory closure is available. Diagnostic values below are again for the REJECTED moderate candidate:

| Quantity | Value |
|---|---:|
| Delta E total, energy subtraction | -5.01578192868557e-9 J/m |
| Delta E, stable BE-work evaluator | -5.015781598995684e-9 J/m |
| Interval trapezoid D_CH | 5.883640541989877e-8 J/m |
| Interval trapezoid D_visc | 1.5503149781893332e-17 J/m |
| Interval trapezoid D_slip | 1.2202676241099858e-29 J/m |
| Interval trapezoid D_total | 5.883640543540192e-8 J/m |
| B_BE (not heat) | 1.3461097811392813e-9 J/m |
| Weak-work defect | -1.1568203584615347e-19 J/m |
| Decomposition roundoff | -3.0810823456375175e-20 J/m |
| Work-identity roundoff | -3.0814338973978525e-20 J/m |
| Continuum local defect (stable Delta E) | 5.3820623836406236e-8 J/m |
| Energy-subtraction budget defect | 5.382062350671635e-8 J/m |

Energy non-growth and every discrete work gate pass. The large continuum quadrature defect is expected to be investigated by the qualified fine clock, not by this dt=1e-6 one-root comparison; it is neither hidden nor used as an algebraic comparison veto. No temporal qualification is claimed here.

The mandatory full-trajectory DeltaE, cumulative total physical dissipation, final budget/|DeltaE|, sum B_BE and max trajectory weak-work values are all NOT RUN.

## 17. Phi coupling

NOT RUN. Primary scale is the initial disturbance norm, not ||phi||. No production full-minus-isolated fields exist.

## 18. Chemical coupling

NOT RUN. The initial/current chemical-disturbance floor and no-mean-subtraction convention remain frozen.

## 19. D_CH coupling

NOT RUN. No scalar full/isolated clock comparison or integrated CH coupling gap exists.

## 20. Kinetic energy

Full trajectory maximum: NOT RUN. Rejected moderate E_kin=2.860256305129656e-17 J/m is diagnostic only, not a transfer measurement.

## 21. Viscous/slip dissipation

Full trajectory instantaneous/cumulative hydro ratios: NOT RUN. Moderate powers in section 8 cannot substitute for accumulated coupling.

## 22. Advective-vs-diffusive phase coupling

At the rejected moderate candidate, independently assembled P2 Riesz norms:
advection 0.0000031828540829583806,
diffusion 6.065134674195071,
ratio 5.247787978229503e-7.
This is diagnostic only; it does not authorize downstream PDE or transfer.

## 23. Transfer relative to isolated temporal uncertainty

| Coupling observable | Measured max/final | Frozen threshold | Pass |
|---|---|---|---|
| phi / initial disturbance | NOT RUN | 1e-3 | NOT RUN |
| F_CH / initial excess | NOT RUN | 1e-3 | NOT RUN |
| Significant D_CH relative | NOT RUN | 1e-2 | NOT RUN |
| E_kin / initial excess | NOT RUN | 1e-4 | NOT RUN |
| Significant hydroD / D_CH | NOT RUN | 1e-4 | NOT RUN |
| Cumulative hydroD / initial excess | NOT RUN | 1e-4 | NOT RUN |
| Max phi coupling | NOT RUN | <=3.691188008762449e-8 | NOT RUN |
| Final integrated-D_CH coupling | NOT RUN | <=7.033312454099817e-11 J/m | NOT RUN |
| Final F_CH coupling | NOT RUN | <=6.0652038946784614e-12 J/m; floor 1e-13 | NOT RUN |
| Full final budget / actual energy-change magnitude | NOT RUN | 5% | NOT RUN |

No threshold, dt, guess, line search or model change followed the STOP.

## 24. Combined error proxy

Full coupling components and combined phi/D uncertainty proxies: NOT RUN. The isolated L0-L2 gaps cannot alone be relabelled full-CHNS uncertainty bounds.

## 25. Computational cost and tests

| Measured moderate attempt category | Wall seconds |
|---|---:|
| Initialization/JIT/initial diagnostics | 11.520139167994785 |
| SNES application timer | 7.387541513002361 |
| Post-step diagnostics | 0.1658894659994985 |
| Certified resolution | 0.08746242199413246 |
| Full BC/continuity/phase audit | 0.08488900800148258 |
| BE work | 0.10541882299730787 |
| Solve plus post-check wrapper | 7.869603218998236 |
| Whole failed scientific session | 19.65485542299575 |

SNES CPU timer: 7.386859311 CPU s, separate from wall time. Rejected archival I/O was not separately instrumented; no invented value is supplied.

Inclusive PETSc event timings (do not sum with application times): residual evaluation 0.48310987 s, Jacobian evaluation 2.443419233 s, numerical LU 3.344331966 s, symbolic LU 1.155882786 s, MatSolve 0.247712619 s. These are process events and include the initial consistent-mu solve where applicable.

Host full pytest: **289 passed, 38 skipped, 2 deselected**. Pinned STEP3–STEP3A8 Docker: **212 passed**, no skipped tests. New tiny tests include arbitrary residual/Jacobian/BC equivalence, strict moderate roots, extreme-dt retained-rate archive, stationary toy, multistep dt change, new-process restart, transactional failures and complete-history coupling/figure integration.

Remote CI for the new source: **NOT RUN; no push**. Existing workflow now includes the new tiny tests; production 96x48 runs are not in CI.

Source remained frozen throughout production and read-only diagnosis:
base HEAD `1cc78a1f4108876a1610570b48138381f29a2484`;
multiphase source `031f988c2aeca347c554f5124eb84ed5c748090670ac0468f458ba9bfc20ae35`;
runner `9d517d735923dcd91431d73de9b260ae9e6e251512f760d25f7896ad76c2bc9d`;
policy `7de1ccc5132b6eaab63d0c7f31a24bfb556fb36b02eee57c1bb7441b65436e25`;
exact source archive `3698902dca348e7d36174308b2c7da7d26da6b78e7ff52fd861adae3886f9b78`.
A later commit is not retrospectively substituted for execution provenance.
Pinned DOLFINx 0.10.0 / PETSc 3.24.0 / MUMPS 5.8.1, float64, one rank.

## 26. Remaining limitations / NOT RUN

Production P1 rate, P2 absolute/rate, exact full target tiny, pilot, production restart, first dt transition, continuous full L0, all trajectory coupling/dominance/closure metrics and combined proxies: NOT RUN.

The seven production coupling PDFs (phi, CH power, CH energy, hydro powers, velocity, budget, solver iterations) were not created without production trajectory data. Their generation was exercised on tiny complete test histories only.

Full CHNS L1/L2, unequal density, body force, BDF2, moving contact line, angle quench, epsilon/mobility/slip sensitivity, falling film, gravity, water-air and tank remain outside scope and were not started. STEP3A6/7 historical evidence and FAIL/PASS records are untouched.

## 27. Scientific verdict

TINY-MESH FULL CHNS PHASE-RATE ALGEBRAIC BRIDGE PASS;

PRODUCTION MODERATE COMPARISON STOPPED BY FROZEN
STRONG-DIVERGENCE NON-REGRESSION GATE;

COUPLING TRANSFER NOT QUALIFIED.

MODEL NOT YET VALIDATED.

The next decision requires a separate assessment of the divergence-envelope applicability at the new moderate interval. This checkpoint does not loosen the gate or authorize a retry.
