# STEP 3A.4 — cost-efficient CH validation checkpoint

**MODEL NOT YET VALIDATED.** Linear architecture succeeds; nonlinear pilot
stops on the FIRST interval. No accepted stiff-bump trajectory, no three-level
qualification and no full-CHNS coupling transfer. All missing metrics below
are explicitly N/A, not linear-reference substitutes.

## 1. Goal

Separate full sparse LINEAR verification, full NONLINEAR isolated CH, and
targeted full NONLINEAR CHNS. Preserve the physical model and scientific gates.
Implementation base: ec01bf3a570cbf2d58c8e16bc133f43f8a7c8baa.
Policies were frozen in STEP3A4_DESIGN.md before new calculations; policy SHA
`f28b662ff5d7d5c55093efffdeeb644d4406b2bb8221921ab41bafcb62e07d6a`. Historical reports/data unchanged.

## 2. Why STEP3A.3 stopped

The reduced-linear passing 1796-step plan implied 4.65 hours level0 and
32.54283029 hours for full-CHNS plan/half/quarter, above the unchanged
3/10-hour old limits. That was a cost checkpoint, not nonlinear failure evidence.

## 3. Old cost estimate

The conservative historical maximum was 9.31865964330 seconds per accepted
full-CHNS interval. It was inappropriately also used to suppress sparse-linear
verification. It is NOT an isolated-CH timing and NOT measured total CPU usage.

## 4. Separation of analysis/isolated/full-CHNS costs

New fixed budgets: analysis 1 h/stage and 2 h total; optimizer 30 min CPU with
finite evaluation/round/split guards; isolated level0 2 h, pilot+three levels
8 h total; targeted CHNS probe 2 h. These roles have separate authorizations.
Actual charged analysis: 278.337140 s. No budget was increased.

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

| stage | sec/accepted step | steps | total wall runtime [s] |
| --- | --- | --- | --- |
| historical full CHNS (projected 3-level series) | 9.31865964 | 12572 | 117154.189 |
| full sparse linear baseline | 0.084830187 | 1796 | 152.355016 |
| full sparse linear optimized | 0.0891180356 | 1068 | 95.178062 |
| isolated failed pilot | undefined: 0 accepted | 0 accepted / 1 attempted | 17.1142374 |
| isolated CH level0 | N/A — NOT RUN | 1068 planned / 0 accepted | N/A — NOT RUN |
| isolated CH level1 | N/A — NOT RUN | 2136 planned / 0 accepted | N/A — NOT RUN |
| isolated CH level2 | N/A — NOT RUN | 4272 planned / 0 accepted | N/A — NOT RUN |
| full CHNS coupling probe | N/A — NOT RUN | 745 prospective / 0 executed | N/A — NOT RUN |

Baseline: 19 numerical setups for 19 distinct dt; setup/factorization
19.463450 s;
solves 66.945816 s;
reference/diagnostics 60.792460 s;
peak RSS 385.641 MiB.
Optimized: 19 setups for
18 distinct dt (adjacent equal-dt blocks
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
exception, so the application failure-safe timer (12.620930 s)
is authoritative for that attempt. See the [PETSc profiling manual](https://petsc.org/release/manual/profiling/).

## 7. Old 898/1796 plan

| plan | steps | endpoint L2 | max L2 | energy error | D2 error | full-sparse pass |
| --- | --- | --- | --- | --- | --- | --- |
| original 898 | 898 | 0.00155535131 | 0.00183649168 | 0.00165580544 | 0.00345324238 | NOT RUN: reduced FAIL |
| uniform 1796 | 1796 | 0.000778338881 | 0.000920446467 | 0.000828515662 | 0.00171727794 | PASS |
| optimized 1068 | 1068 | 0.000787203393 | 0.000948259981 | 0.000854204556 | 0.00216351617 | PASS |

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

| candidate | action | block | steps | endpoint | max L2 | energy | D2 | PASS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | original | — | 898 | 0.00155535131 | 0.00183649168 | 0.00165580544 | 0.00345324238 | False |
| 1 | trial_split | 17 | 1032 | 0.000835681308 | 0.0011555455 | 0.00103905098 | 0.00266428099 | False |
| 2 | trial_split | 18 | 899 | 0.00155520891 | 0.00183649168 | 0.00165580544 | 0.0034532423 | False |
| 3 | accept_split | 17 | 1032 | 0.000835681308 | 0.0011555455 | 0.00103905098 | 0.00266428099 | False |
| 4 | trial_split | 16 | 1068 | 0.000787203393 | 0.000948259981 | 0.000854204556 | 0.00216351617 | True |
| 5 | trial_split | 17 | 1300 | 0.000475192644 | 0.00108872217 | 0.000960167296 | 0.00227239595 | False |
| 6 | trial_split | 18 | 1033 | 0.000835539636 | 0.0011555455 | 0.00103905098 | 0.00266428091 | False |
| 7 | accept_split | 16 | 1068 | 0.000787203393 | 0.000948259981 | 0.000854204556 | 0.00216351617 | True |
| 8 | trial_merge | 16 | 1032 | 0.000835681308 | 0.0011555455 | 0.00103905098 | 0.00266428099 | False |

All nine tested/accepted candidate records are retained. Only blocks 17 and
16 were refined, in that order. The result is a legal one-block merge fixed
point. Independent repeat semantic hashes match. Optimization and repeat took
4.158055 s including setup, well inside guards.
The repeat compares the deterministic optimizer CORE semantic SHA
`2c22ab5db535cefe83a7043cb45da463e9385c3d80de34bbed2a3cd45f358d45` (fresh search history/cache per invocation).
The frozen execution-artifact SHA below additionally includes archived source
provenance. Gzip archive timestamps are not part of optimizer reproducibility.
All six scientific source archives were independently checked against both
recorded source hashes and archive-file hashes; all match.

## 10. Optimized plan

1068 intervals: 728 fewer, **40.534521% shorter** than 1796. Half/quarter
would be 2136/4272 intervals with identical boundaries. dt0=4.657657028534679e-12 s.
Plan semantic SHA: `5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8`.

| block (0-based) | t_start | t_end | dt | steps |
| --- | --- | --- | --- | --- |
| 0 | 0 | 5.68234157e-10 | 4.65765703e-12 | 122 |
| 1 | 5.68234157e-10 | 1.08057643e-09 | 9.31531406e-12 | 55 |
| 2 | 1.08057643e-09 | 1.93758532e-09 | 1.86306281e-11 | 46 |
| 3 | 1.93758532e-09 | 3.3162518e-09 | 3.72612562e-11 | 37 |
| 4 | 3.3162518e-09 | 7.34046748e-09 | 7.45225125e-11 | 54 |
| 5 | 7.34046748e-09 | 1.44946287e-08 | 1.49045025e-10 | 48 |
| 6 | 1.44946287e-08 | 2.96972212e-08 | 2.9809005e-10 | 51 |
| 7 | 2.96972212e-08 | 4.93711645e-08 | 5.961801e-10 | 33 |
| 8 | 4.93711645e-08 | 1.05412094e-07 | 1.1923602e-09 | 47 |
| 9 | 1.05412094e-07 | 2.24648114e-07 | 2.3847204e-09 | 50 |
| 10 | 2.24648114e-07 | 4.10656305e-07 | 4.7694408e-09 | 39 |
| 11 | 4.10656305e-07 | 8.49444858e-07 | 9.53888159e-09 | 46 |
| 12 | 8.49444858e-07 | 1.1928446e-06 | 1.90777632e-08 | 18 |
| 13 | 1.1928446e-06 | 2.75722118e-06 | 3.81555264e-08 | 41 |
| 14 | 2.75722118e-06 | 4.51237539e-06 | 7.63110528e-08 | 23 |
| 15 | 4.51237539e-06 | 7.10695118e-06 | 1.52622106e-07 | 17 |
| 16 | 7.10695118e-06 | 1.80957428e-05 | 1.52622106e-07 | 72 |
| 17 | 1.80957428e-05 | 9.99011913e-05 | 3.05244211e-07 | 268 |
| 18 | 9.99011913e-05 | 0.0001 | 9.88086651e-08 | 1 |

## 11. Reduced reference verification

The qualified historical rational M-Krylov reference, matrices and perturbation
are reused, not regenerated or relabelled as newly generated evidence. Inputs
and reference arrays are fingerprint checked. Reduced optimized gates PASS.

## 12. Full-sparse optimized verification

| metric | reduced optimized | full sparse optimized | full minus reduced |
| --- | --- | --- | --- |
| endpoint_relative_L2_error | 0.000787203393 | 0.000787203392 | -4.21131879e-13 |
| max_checkpoint_relative_L2_error | 0.000948259981 | 0.000948259981 | 8.16432425e-14 |
| max_quadratic_energy_relative_error | 0.000854204556 | 0.000854204556 | -1.87447171e-13 |
| integrated_D2_relative_error | 0.00216351617 | 0.00216351642 | 2.49412713e-10 |
| max mass/domain | not independently assessed here | 1.90977954e-18 | N/A |

All five full-sparse gates PASS. Actual runtime 95.178062 s,
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
dt=4.657657028534679e-12 s. Accepted intervals=0. SNES reason=-6
(DIVERGED_LINE_SEARCH), 11 iterations; final residual=1.11091103673e-08.
Initial residual=0.784388568715275. With unchanged atol=1e-10, rtol=1e-9,
stol=0, the residual target is max(atol, rtol*initial)=7.84388568715e-10.
The stored reassembly reproduces the reported residual. No retry, smaller dt,
tolerance change or solution correction was attempted. See the
[SNES tolerance definitions](https://petsc.org/release/manualpages/SNES/SNESSetTolerances/).

At t=0: phase mass=-0.6060961027499487, mass/target-domain error
5.51719164e-13, certified cells
8.32813107583, two wall crossings.
D_CH(0)=0.11400313905309394 W/m and excess
5.5428884676578249e-08 J/m reproduce historical values exactly.
Initial phi SHA and psi SHA are unchanged. The DOF range is
[-0.979080699053, 1.01569533975]; no clipping.

Initial phi SHA: `6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681`.
Psi SHA: `831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0`.
Prepared equilibrium: `df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc`.
Mesh: `93dbd2cf76148bcd50afef9e353c55f92c15789ad0e303dbb97435a4672b9b15`. Quadrature degree 12, P2;
all physical/config parameters are loaded unchanged from the prepared cache.

Pilot wall time 17.114237 s, of which failed
SNES attempt 12.620930 s.
This is NOT seconds per accepted interval. Isolated throughput and full-series
cost cannot be established from zero accepted steps. Cost authorization remains false.

## 15. Half refinement

NOT RUN, blocked by pilot failure. No alternate startup or dt retry.

## 16. Quarter refinement

NOT RUN. No formal three-level convergence evidence exists.

| level | min dt (planned) | steps planned/accepted | integral D_CH | DeltaF | final defect | defect/DeltaF | endpoint difference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 4.65765703e-12 | 1068/0 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |
| 1 | 2.32882851e-12 | 2136/0 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |
| 2 | 1.16441426e-12 | 4272/0 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |

## 17. Integrated D_CH convergence

Not measured: no accepted isolated intervals. No Richardson extrapolation,
linear D2 integral or rejected-state endpoint power substitutes for the missing
nonlinear physical-dissipation integrals. Last relative D_CH difference: N/A.

## 18. Energy closure

Nonlinear global closure and successive defect decrease: N/A. The initial
energy is 0.065642782684546278 J/m. A near-zero work residual for one
REJECTED Newton iterate (below) does not qualify a nonlinear time history.

## 19. Endpoint phi convergence

N/A: no level0/half/quarter endpoints. No conditional transfer is inferred.

## 20. Early BE work decomposition

There are NO accepted BE-work intervals. The following is a read-only
postmortem of the stored REJECTED iterate, clearly not trajectory evidence:

| n | dt | DeltaF | Dold | Dnew | B_BE | weak defect | trapezoid gap | continuum defect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 REJECTED | 4.65765703e-12 | -5.23519646e-13 | 0.114003139 | 0.111872687 | 2.45501711e-15 | -2.7945577e-20 | 4.96145805e-15 | 2.50641805e-15 |

Decomposition roundoff=5.05740801e-21 J/m.
Physical continuum dissipation is D_CH+D_visc+D_slip. B_BE is the signed
NUMERICAL BE work remainder, diagnostic-only, never physical heat.

Residual localization (algebraic norms are solver-conditioning diagnostics,
not scientific stationarity norms): phase=1.11091103673e-08;
chemical=5.84182476984e-17. The time and flux vector norms
are ~0.773630519 and nearly cancel. Evaluating the same time term using
M0*(phi_coeff-old_coeff)/dt differs from the literal UFL evaluation by
8.91767199985e-09. Its complete phase residual
is still 1.19088969405e-08; merely reassembling
the rejected coefficients is NOT a demonstrated fix. Stable-time UFL vs sparse
assembly differs by 4.12945446052e-16.

The L2 discrepancy between subtracting interpolated fields and interpolating
coefficient differences is 3.64692778999e-17;
division by tiny dt amplifies it to 7.82996207675e-06.
Phase Riesz L2 residual=2.83642501247e-06.
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
prospective prefix estimate is 1.928445 h,
below the 2 h probe limit; this does not bypass the scientific prerequisite.
Coupling thresholds remain phi/E 1e-3, D_CH 1e-2 above 1e-8*D_CH(0),
kinetic/excess and hydroD/D_CH 1e-4, mass 1e-10, certified cells 8.

## 23. Full CHNS coupling results

| prospective time | phi L2 rel | H1 rel | energy rel | D_CH rel | E_kin/excess | hydroD/D_CH |
| --- | --- | --- | --- | --- | --- | --- |
| 4.65765703e-10 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |
| 4.86967358e-07 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |
| 2.45197697e-06 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |
| 9.85414908e-06 | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN | N/A — NOT RUN |

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
full-sparse runtime reduced from 152.355016 to
95.178062 s. Hypothetical full-CHNS plan/half/quarter
would drop from 32.542830 to 19.351750 wall hours, still not authorized.
New measured scientific stages consumed 295.451378 wall s and 282.194768 process
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

Local pure suite: 263 passed, 18 skipped
(FEniCSx dependencies), 2 validation-marked deselected. Pinned Docker workflow:
166 passed, 0 failures,
0 errors. The initial Docker invocation
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
