# STEP 3A.7 measured report

## 1. Goal and inherited evidence

STEP3A6 proved one tiny BE step and absolute-state cancellation, not multi-step temporal convergence. This iteration tests nonlinear isolated CH only. Historical STEP3A2–6 files and FAIL/checkpoint verdicts are unchanged.

## 2. Immutable physics/source

Execution base: `dfe272049123ecbd9d200af8ff070ae4eb4b08d3`; source `e762d8a21f3b0b35c72a49efa5cff581f0153d9bebc40c0c266ecb31265808f9`; runner `6fef8605de514ba61de71a49e129a84177c803b5dd7f87f17627ebc0f8e64e27`; policy `1113da00ba12ddc746e362841eeee44c20a6ee87422342fe6e78b9c798dc2a47`. These identify the dirty execution snapshot, not the later result commit. Inputs, P2 mesh, quad12, physical parameters and original production Newton controls remain fixed.

| input | SHA |
| --- | --- |
| equilibrium_fingerprint | df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc |
| initial_phi_sha256 | 6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681 |
| optimized_plan_sha256 | 5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8 |
| perturbation_coefficient_sha256 | 831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0 |


## 3. Level0 schedule and derived levels

| level | intervals | derived schedule SHA |
| --- | --- | --- |
| 0 | 1068 | 21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31 |
| 1 | 2136 | a289d4ea7769f6d0be6e06d5bd4b1aad5c587696c9303f7d05e68a68e9753a08 |
| 2 | 4272 | 9c128247e195c51d87dacec7f93469103a95389df2f348577f4ae142d5b0ab72 |

Parent plan hash is unchanged. Each parent endpoint is retained exactly; child solver dt=parent dt/factor. Clock-rounding discrepancy is recorded, not used to adapt dt.

## 4. Multi-step phase-rate runner

Unknowns are retained (r_phi,mu); phi_new=phi_old+dt*r_phi. First guess is semidiscrete rate; every subsequent guess is the previous accepted rate, unscaled. All intervals use atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS. No retry or absolute/BDF2 switch.

## 5. Transaction/checkpoint design

Private physical candidate and copied persistent Diagnostics pass material/finite/mass/certified-resolution/topology/energy/work checks before publication. True rate, state/old/older, time/index, maps, diagnostic initial/previous/cumulative and discrete work sums are archived with hashes. Journal prefixes and COMPLETE markers are checked; rejected states never restart. Scalar and timing journals are append-only; selected fields/checkpoints are atomically published.

## 6. Pilot

| metric | result |
| --- | --- |
| accepted intervals | 127 |
| final pilot time | 6.14810727767e-10 |
| first block boundary index | 122 |
| crossed boundary? | True |
| max Newton | 1 |
| max mass/domain | 1.97372982156e-14 |
| min certified cells | 8.32813107583 |
| max weak work | 8.96446929553e-20 |
| wall time including fork/init | 235.081584964 |
| sec/accepted interval | 1.58242484494 |
| projected L0 wall s | 3159.66526613 |

Pilot is a continuous reusable L0 prefix, not temporal qualification. Its
`execution_status=incomplete` refers to the then-unfinished 1068-step L0,
not to a failed pilot: all 127 required prefix intervals passed and L0 was
authorized. The same dataset was subsequently continued to completion.

## 7. First block transition

The primary transition passed. The five preceding Newton counts were all 1;
the first post-boundary interval also required 1 iteration. The previous
accepted rate was retained without rescaling, including its coefficient SHA.

| first post-boundary quantity | value |
| --- | --- |
| step / physical endpoint [s] | 123 / 5.775494715383002e-10 |
| old dt / new dt [s] | 4.657657028534679e-12 / 9.315314057069358e-12 |
| SNES reason / KSP reason | 2 / 4 |
| initial / final residual | 4.15548491130545e-9 / 1.0118177987671685e-14 |
| effective target | 1e-10 |
| guess rate L2 / converged rate L2 | 102.1389234485376 / 101.39968724130613 |
| dt max(abs(rate)) | 2.895378395728337e-8 |
| dt reductions | 0 |

Full controls, trace and hashes are in
[block_transition/status.json](validation_results/step3a7/block_transition/status.json).
The optional alternate-guess boundary control was not elected in the design.

## 8. Restart equivalence

PASS. A fresh diagnostic process resumed the accepted production checkpoint
at step 25 and independently executed steps 26–28. At all three endpoints,
the measured phi, mu and retained-rate relative L2 differences, cumulative-D
relative difference, and absolute F/budget differences were exactly 0 in the
stored comparison. This fork did not replace primary fields. The separate
tiny-process regression additionally restarted across a dt change.

See [pilot/restart_check.json](validation_results/step3a7/pilot/restart_check.json).

## 9. Cost qualification

Committed machine limits: 7200s L0, 28800s whole isolated series. Forecast=1.25*max(mean,p95) of measured accepted phase-rate intervals, plus measured spent/setup. No historical CHNS seconds/step. Pilot reuse is charged once. {"authorized": true, "level0_cost_ok": true, "mean_s": 1.5824248449448926, "p95_s": 2.474090010199916, "projected_level0_wall_s": 3159.6652661286516, "projected_s_per_step": 3.092612512749895, "projected_series_wall_s": 23005.996861163978, "remaining_counts": [941, 2136, 4272], "remaining_wall_s": [2924.5836811646514, 6620.255633900776, 13226.075961134553], "safety_factor": 1.25, "series_cost_ok": true, "spent_wall_s": 235.08158496400029}

## 10. Level0

| metric | L0 | L1 | L2 |
| --- | --- | --- | --- |
| execution_status | complete | complete | complete |
| accepted_steps | 1068 | 2136 | 4272 |
| min_dt | 4.65765702853e-12 | 2.32882851427e-12 | 1.16441425713e-12 |
| max_dt | 3.05244211022e-07 | 1.52622105511e-07 | 7.63110527555e-08 |
| wall_s | 1414.1444098 | 2646.92603727 | 5418.65345229 |
| max_Newton | 1 | 1 | 1 |
| max_mass_domain | 4.63672310423e-13 | 4.59354776439e-13 | 4.45168593346e-13 |
| max_target_mass_domain | 6.34369100459e-13 | 6.43620958998e-13 | 6.38069843875e-13 |
| min_certified_cells | 8.32813107583 | 8.32813107583 | 8.32813107583 |
| final_F | 0.0656427394272 | 0.0656427394232 | 0.0656427394212 |
| integrated_D_CH | 4.33589839823e-08 | 4.33120093054e-08 | 4.32886508578e-08 |
| DeltaF | -4.32573141196e-08 | -4.32613577322e-08 | -4.32633793235e-08 |
| final_budget_defect | 1.01669862723e-10 | 5.06515731348e-11 | 2.52715342869e-11 |
| budget_over_DeltaF | 0.00235035079713 | 0.0011708271721 | 0.000584132231049 |

Only COMPLETE levels enter formal convergence. No partial final horizon is silently substituted.

## 11. Level1

Status: complete. Fresh identical t=0; every L0 interval split exactly in two. See level table.

## 12. Level2

Status: complete. Fresh identical t=0; every L0 interval split exactly in four. No extrapolation substitutes for this level.

## 13. Mass and resolution over all levels

Both initial and target mass/domain are checked every accepted state, as are certified P2 cells>=8 and two wall crossings. Full minima/worst-cell witnesses are in summary.json. Figure: analysis/resolution.pdf.

## 14. D_CH trajectories

Independent physical power includes t=0 and every interval. Figures: dissipation.pdf and dissipation_log_time.pdf. Isolated CH has u=0 by construction; no measured full-CHNS zero powers are claimed.

## 15. Integrated D_CH convergence

| observable | L0–L1 | L1–L2 | passes? |
| --- | --- | --- | --- |
| final phi L2 | 2.0667419446e-08 | 1.03399527806e-08 | True |
| max common-time phi L2 | 2.45892706579e-08 | 1.23226176604e-08 | diagnostic |
| integrated D relative | 0.00108456471252 | 0.00053959749628 | True |
| final F difference | 4.04361266693e-12 | 2.02159122775e-12 | True |
| absolute final budget [J/m] | 1.01669862723e-10 → 5.06515731348e-11 | 5.06515731348e-11 → 2.52715342869e-11 | True |

Physical integral uses trapezoids, never -DeltaF. The successive relative
differences are 0.108456471252% and 0.053959749628%; both decrease and the
last is below the unchanged 5% threshold.

| interval contribution / total physical I_D | L0 | L1 | L2 |
| --- | --- | --- | --- |
| first interval (fraction, not percent) | 1.21318816892e-5 | 6.10071243358e-6 | 3.05921975786e-6 |
| largest interval (fraction, not percent) | 0.00526480738240 | 0.00264361903968 | 0.00132465657767 |

These fractions are diagnostics, not added acceptance thresholds.

## 16. Energy-budget convergence

Budget=F(t)-F(0)+integrated D_CH. Shared positive-scale/change/symmetric normalizations are retained; finest nonzero abs(B)/abs(DeltaF)<=5% is separately required without hiding it behind a significance floor. Figure: energy.pdf.

{"nonzero_relaxation": true, "dissipation_gaps_decrease": true, "last_D_difference": true, "budget_decreases": true, "finest_closure": true, "endpoint_phi_converges": true, "endpoint_energy_converges": true, "physical_all_levels": true, "shared_finest_energy_gate": true}

Finest final closure is **0.058413223105%**, below the unchanged 5% limit.
The shared energy-change applicability flag is false because its inherited
significance floor is 6.56427826845e-8 J/m. Its legacy numeric zero placeholder
therefore does NOT mean zero defect. The explicitly required actual-change
ratio above is measured independently of that applicability flag.

## 17. Endpoint phi convergence

Consistent P2 mass norm gaps: [2.0667419445980026e-08, 1.0339952780595003e-08]. Both field and dissipation convergence are mandatory, not energy-only evidence.

## 18. Common-time convergence

Max phi gaps: [2.4589270657912066e-08, 1.2322617660374471e-08]. Exact shared accepted endpoints; no interpolation. Phi/mu/F/D/integral comparisons are in convergence.json and phi_convergence.pdf.

## 19. Observed temporal orders

{"phi": 0.999128666387107, "D_gaps": 1.0079388419241038, "budget": [1.0052131126823791, 1.0030938209656557]}. Diagnostic only; no demand for exactly first order and no fitted thresholds.

## 20. BE work diagnostics

| interval | dt | Delta_E | Dold | Dnew | B_BE | weak defect | trap gap | continuum defect |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L0:1 | 4.65765702853e-12 | -5.23519657765e-13 | 0.114003139053 | 0.11187268676 | 2.45501710735e-15 | -3.94251501092e-20 | 4.96145804793e-15 | 2.50640607228e-15 |
| L0:2 | 4.65765702853e-12 | -5.13898748375e-13 | 0.11187268676 | 0.109828345479 | 2.35602831076e-15 | 4.93046236595e-20 | 4.76092026781e-15 | 2.40493714684e-15 |
| L0:3 | 4.65765702853e-12 | -5.04664038937e-13 | 0.109828345479 | 0.107865842565 | 2.26194259308e-15 | 1.34999975521e-21 | 4.5703327458e-15 | 2.3083935709e-15 |
| L0:4 | 4.65765702853e-12 | -4.95796261604e-13 | 0.107865842565 | 0.105981145214 | 2.17248510523e-15 | 4.62587175531e-20 | 4.3891369329e-15 | 2.21670122542e-15 |
| L0:5 | 4.65765702853e-12 | -4.8727765129e-13 | 0.105981145214 | 0.104170445764 | 2.0873978207e-15 | -4.66174198011e-20 | 4.21680850859e-15 | 2.12936609822e-15 |
| L1:1 | 2.32882851427e-12 | -2.63600917643e-13 | 0.114003139053 | 0.112921352874 | 6.2652668776e-16 | 7.19239105283e-20 | 1.25964725021e-15 | 6.3319604917e-16 |
| L1:2 | 2.32882851427e-12 | -2.61120908313e-13 | 0.112921352874 | 0.111861927872 | 6.13591998006e-16 | -6.78470363241e-20 | 1.23360957686e-15 | 6.19948552059e-16 |
| L1:3 | 2.32882851427e-12 | -2.58691748354e-13 | 0.111861927872 | 0.110824294327 | 6.00986566477e-16 | 9.03651425554e-21 | 1.20823529299e-15 | 6.07263641685e-16 |
| L1:4 | 2.32882851427e-12 | -2.56312477941e-13 | 0.110824294327 | 0.109807898992 | 5.88700881292e-16 | -7.75694041736e-21 | 1.18350521938e-15 | 5.94793541876e-16 |
| L1:5 | 2.32882851427e-12 | -2.53981709349e-13 | 0.109807898992 | 0.108812204567 | 5.76725730792e-16 | -1.53082987856e-20 | 1.15940078444e-15 | 5.82656130539e-16 |
| L2:1 | 1.16441425713e-12 | -1.32270381831e-13 | 0.114003139053 | 0.113457998912 | 1.58274764481e-16 | 1.9661182723e-21 | 3.17384476103e-16 | 1.59114164171e-16 |
| L2:2 | 1.16441425713e-12 | -1.31640554128e-13 | 0.113457998912 | 0.112918558059 | 1.56621077028e-16 | 4.71587915386e-20 | 3.14066310245e-16 | 1.57491080718e-16 |
| L2:3 | 1.16441425713e-12 | -1.31017403568e-13 | 0.112918558059 | 0.112384742897 | 1.5498873515e-16 | -1.99272573274e-20 | 3.10790992738e-16 | 1.55784338151e-16 |
| L2:4 | 1.16441425713e-12 | -1.3040068335e-13 | 0.112384742897 | 0.111856480908 | 1.53377425092e-16 | -2.41943906456e-20 | 3.07557895566e-16 | 1.54155667662e-16 |
| L2:5 | 1.16441425713e-12 | -1.29790351379e-13 | 0.111856480908 | 0.111333700639 | 1.51786838953e-16 | -1.75496675914e-20 | 3.043663994e-16 | 1.52563343847e-16 |

Signed numerical B_BE is a discrete BE work remainder, NOT physical heat. It is never added to physical dissipation. Basic work is measured every step; expanded same-state residual audit is sparse and read-only, never an absolute-residual acceptance veto.

| whole-trajectory work diagnostic [J/m] | L0 | L1 | L2 |
| --- | --- | --- | --- |
| max absolute weak-work defect | 6.00161720618e-19 | 7.78117632059e-19 | 6.87090257207e-19 |
| max absolute decomposition roundoff | 4.61362352426e-19 | 5.44111586770e-19 | 4.61415176689e-19 |
| sum dt D_new | 4.31566728413e-8 | 4.32109296984e-8 | 4.32381300449e-8 |
| sum B_BE | 1.00616603979e-10 | 5.04048246917e-11 | 2.52266319395e-11 |
| sum weak-work defect | 2.14469912526e-19 | 5.95054734691e-19 | -4.23530807444e-19 |

All absolute work limits remain 1e-12 J/m. The summed stable interval Delta_E
differs from total F subtraction by 2.47e-14 / 2.32e-14 / 2.26e-14 J/m;
both evaluations are retained, and the primary continuum budget still uses F.

## 21. Solver statistics

Newton distributions, maximum/index, KSP failures, rejected count, all dt transitions and rate/guess norms are in summary.json. Figures: solver_iterations.pdf and phase_rate.pdf.

All 7476 primary accepted intervals required exactly **one Newton iteration**:
min/mean/p50/p95/max are all 1 on every level. Thus the maximum is tied at
every step (the summary records its first occurrence, step 1). KSP failures
and rejected intervals are 0. All finite/material/topology checks passed,
with two wall crossings; no positive free-energy increment was measured.

First actual dt-change indices (before / first-after / second-after Newton
counts are 1 / 1 / 1 wherever the second interval exists):

| level | dt-change indices |
| --- | --- |
| L0 | 123, 178, 224, 261, 315, 363, 414, 447, 494, 544, 583, 629, 647, 688, 711, 800, 1068 |
| L1 | 245, 355, 447, 521, 629, 725, 827, 893, 987, 1087, 1165, 1257, 1293, 1375, 1421, 1599, 2135 |
| L2 | 489, 709, 893, 1041, 1257, 1449, 1653, 1785, 1973, 2173, 2329, 2513, 2585, 2749, 2841, 3197, 4269 |

Adjacent original blocks with identical dt are not falsely labelled dt changes.

Local test records:

| run | tests | failed | errors | skipped |
| --- | --- | --- | --- | --- |
| docker_final.xml | 199 | 0 | 0 | 0 |
| docker_preflight.xml | 199 | 0 | 0 | 0 |
| preflight_trajectory.xml | 6 | 1 | 0 | 0 |
| preflight_trajectory_r2.xml | 6 | 0 | 0 | 0 |
| pure_final.xml | 314 | 0 | 0 | 30 |
| pure_preflight.xml | 314 | 0 | 0 | 30 |

The early `preflight_trajectory.xml` failure was a test-harness subprocess
PYTHONPATH issue, fixed and retested BEFORE source freeze; it is retained as
development history. Both full preflight suites and both final suites pass.
Final host result: 284 passed, 30 expected container-only skips, 2 deselected
(41.82s). Final pinned Docker result: 199 passed (53.37s). The Docker warning
is the existing deliberately under-resolved negative test, not this production
mesh. No implementation change occurred after the scientific pilot.

Remote GitHub CI NOT RUN (no push). Production 96x48 is excluded from CI.

## 22. Computational cost

Measured scientific wall time including initialization/fork: 9492.55s. Not CPU time. Stage/session records separate initialization/JIT, inclusive PETSc events, SNES, diagnostics, resolution, work, sparse residual audit and archival I/O. Figure: cost.pdf. Timing-record write itself is excluded from per-step sample but included in session wall time.

Total: **2 h 38 min 12.55 s**, below the frozen 8 h limit; local tests and
report generation are separate from scientific trajectory wall time.
L0 forecast was 3159.665 s (52.66 min), compared with 1414.144 s primary
trajectory wall time (23.57 min), or 1426.971 s including its diagnostic fork.
The full L0/1/2 accepted-step means are 1.29032 / 1.22893 / 1.25936 s.
After L0 the remaining-series forecast passed at 15336.782 s total, and
after L1 it passed at 13178.751 s total. Budgets were never increased.

## 23. Remaining limitations

NOT RUN / not qualified: full-CHNS rate solver, CHNS algebraic/multi-step verification, targeted coupling probe, moving 90→60 contact line, theta90/120, epsilon/slip/mobility sensitivity, falling-film PDE, gravity/tank/water-air. Isolated results do not establish full CHNS transfer.



## 24. Scientific verdict

NONLINEAR ISOLATED CH STIFF TRANSIENT TEMPORALLY QUALIFIED ON THE CERTIFIED MESH; FULL CHNS TRANSFER AND REMAINING STEP 3A PHYSICAL BENCHMARKS STILL REQUIRED.

**MODEL NOT YET VALIDATED**.
