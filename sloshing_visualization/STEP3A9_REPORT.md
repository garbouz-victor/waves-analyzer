# STEP3A9 — divergence policy audit and full-CHNS transfer

Overall: **MODEL NOT YET VALIDATED**.

Audit, moderate/tiny bridge, pilot, restart and first dt transition passed.
**STOP: the frozen full-L0 cost forecast exceeds 14400s.** The accepted
continuous prefix has 127/1068 intervals; no full COMPLETE marker exists.
STEP3A8 remains rejected under its own policy. No threshold or source was
changed to continue past the cost gate.

## 1. Goal and immutable inputs

The missing qualification is isolated-to-full hydrodynamic coupling, not a
repeat of isolated temporal refinement. STEP3A7 L0/L1/L2 COMPLETE markers,
history chains and field hashes were verified read-only. L0 has 1068 intervals,
SHA `21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31`.
Initial phi, psi and prepared equilibrium retain their historical hashes.
M0 initial D_CH=0.11400313905309394 W/m and certified cells=8.328131075834218.
No solver core, physical parameters, FE family, time schedule or quadrature changed.

## 2. Mathematical meaning of continuity

The continuous P2-vector / continuous P1 Taylor–Hood pair enforces
`(q_h, div u_h)=0` for pressure-space test functions, to measured nonlinear
tolerance. The divergence is cellwise linear but generally discontinuous:
it need not be in continuous P1. Thus weakly divergence-free is **not**
pointwise divergence-free. This distinction also underlies the reconstruction
of discretely into exactly divergence-free fields in
[Lederer, Linke, Merdon and Schöberl](https://arxiv.org/abs/1609.03701).
No such reconstruction or FE modification is applied here.

Independently assembled `M_Q c=b`, with `b_i=(q_i,div u)`, defines the L2
projection. A separate pressure Function is populated with c and its squared
norm is integrated with quad12. The existing observer's sparse Riesz norm is
computed separately. The complement is integrated as `(div u-P_Q div u)^2`.

## 3. Read-only STEP3A8 rejected-candidate audit

No Newton solve or accepted-state publication took place in this experiment.
The archived residual Function and PETSc iterate agree coefficientwise.

| quantity | measured |
| --- | --- |
| strong divergence L2 | 3.814612027895975e-6 |
| independent projected divergence L2 | 1.1783432940660683e-21 |
| existing continuity Riesz | 1.1783432940661996e-21 |
| projection/Riesz absolute difference | 1.3127920909983307e-34 |
| orthogonal divergence L2 | 3.8146120278959753e-6 |
| projection fraction | 3.0890252677046345e-16 |
| Pythagorean relative squared-norm error | 2.2205452823528433e-16 |

Almost all measured strong divergence is in the pressure-space orthogonal
complement. This conclusion is numerical evidence, not an assumption.

## 4. Cellwise localization

For 18432 cells, sum of cell integrals of divergence squared is
1.4551264923382768e-11; maximum cell contribution is 1.2825469058445928e-13.

| percentile | cell contribution |
| --- | --- |
| p50 | 7.587429430153372e-27 |
| p90 | 1.5689278282312555e-17 |
| p95 | 1.384416162531867e-15 |
| p99 | 2.068374937732617e-14 |

Top 1%, 5%, 10% cells contribute respectively 59.2561%, 97.7524%, 99.9832%
of the **squared** norm. Corresponding norm shares are 76.9780%, 98.8698%,
99.9916%. These two meanings of fraction are not interchangeable.

The predeclared Bernstein-range interface set contributes 4.13648e-9 of
the squared norm; wall-near contributes 9.05133e-8; contact-near 8.07655e-14;
bulk contributes 0.9999999053502688. Interface/wall/contact sets overlap;
bulk excludes both interface and wall-near sets. The map localizes the dominant
contributions around the two bulk bumps, not the diffuse contact region.

Figures: [cell contributions](validation_results/step3a9/divergence_audit/rejected_step3a8/strong_divergence_cells.pdf),
[pressure projection](validation_results/step3a9/divergence_audit/rejected_step3a8/divergence_projection.pdf).
The plotted vertex contour is visual only, never a certified-resolution gate.

## 5. Tiny discrete-kernel demonstration

A deterministic vector is projected into ker(B) of the actual tiny P1/P2
divergence matrix, with the original constrained velocity DOFs excluded.
Kernel dimension=95, coefficient norm=1, algebraic ||Bu||=2.64682e-17,
continuity Riesz=2.72799e-16, strong divergence=1.6561576147827053.
The original state is restored exactly afterwards. This is mathematical FE
evidence, not a production trajectory or a correction of any production state.

## 6. Same-physics spatial study

All initial states passed preparation and resolution **before** transient solves.
M0 loads exact historical coefficients. M1/M2 solve their own stationary
mass-constrained equilibrium, from the same analytic cap seed, with unchanged
equilibrium algorithm and policy. The common target mass is
-0.6060961027495515. The same analytic two-bump definition and deterministic
zero-discrete-mass construction precede each solve, with delta=1e-3.
No solved velocity is transferred between meshes. The larger allocation guard
on refined meshes is metadata, not a physical or nonlinear-accuracy change.

Each row below is one fresh **absolute** BE step on [0,1e-6] with P1 controls:
atol=rtol=1e-12, stol0, max_it30, newtonls/bt, preonly/LU/MUMPS.

| metric | M0 | M1 | M2 |
| --- | --- | --- | --- |
| nx x nz | 96 x 48 | 120 x 60 | 144 x 72 |
| h=max geometric cell diameter [m] | 0.012500000000000178 | 0.010000000000000009 | 0.008333333333333526 |
| crossed triangular cells | 18432 | 28800 | 41472 |
| velocity DOFs | 74306 | 115922 | 166754 |
| pressure DOFs | 9361 | 14581 | 20953 |
| strong divergence L2 | 3.814612027895975e-6 | 2.9097722408663535e-6 | 2.344110670350525e-6 |
| projected divergence L2 | 1.1783432940660683e-21 | 3.774122517275919e-22 | 4.551886037401198e-22 |
| orthogonal divergence L2 | 3.8146120278959753e-6 | 2.9097722408663535e-6 | 2.344110670350525e-6 |
| weak continuity Riesz | 1.1783432940661996e-21 | 3.7741225172788624e-22 | 4.551886037405062e-22 |
| projection fraction | 3.0890252677046345e-16 | 1.2970508358936762e-16 | 1.9418392207226866e-16 |
| grad-u L2 | 4.05648057316023e-6 | 3.141757111679265e-6 | 2.557344300832734e-6 |
| symgrad-u L2 | 3.93740393938612e-6 | 3.027987140205364e-6 | 2.453045526263989e-6 |
| velocity L2 | 7.563407043299337e-9 | 4.523156130722014e-9 | 3.035114530772128e-9 |
| E_kin [J/m] | 2.860256305129656e-17 | 1.0229470691438611e-17 | 4.605960107444485e-18 |
| D_visc [W/m] | 3.100629956378667e-11 | 1.833741224249811e-11 | 1.2034864707847541e-11 |
| certified cells, initial and final | 8.328131075834218 | 10.410163844792772 | 12.492196613751299 |
| SNES reason / Newton iterations | 2 / 2 | 2 / 2 | 2 / 2 |

M0 reproduces the historical candidate's divergence exactly. On refined meshes
the consistent initial D_CH is 0.13787096762514287 and 0.11285495583578013 W/m;
these mesh-dependent discrete values were not forced to match M0.
All three final states pass every non-strong-divergence physical gate.

## 7. Predeclared decision

| criterion | measured | required | pass |
| --- | --- | --- | --- |
| weak continuity, all meshes | max 1.17835e-21 | <=1e-12 | PASS |
| independent projection / Riesz | assembly-roundoff agreement on all | frozen relative + assembly tolerance | PASS |
| Pythagorean consistency | all roundoff-consistent | <=1e-11 relative squared norm | PASS |
| strong-divergence ordering | D_M2 < D_M1 < D_M0 | strict decrease | PASS |
| measured log-log slope | 1.2013838200148381 | >=0.5 | PASS |
| non-divergence physical gates | all pass | all pass | PASS |
| SNES/KSP | all positive; no failure/retry | success | PASS |

Pairwise divergence ratios: 0.7627963786585386, 0.8055993652797345.
This is a measured refinement slope on three close meshes, not an optimal
theoretical convergence-order claim. [Spatial convergence figure](validation_results/step3a9/divergence_audit/divergence_spatial_convergence.pdf).

## 8. Old and new policies

| quantity | historical STEP3A8 | new STEP3A9 |
| --- | --- | --- |
| weak continuity Riesz | hard <=1e-12 | hard <=1e-12 |
| strong divergence L2 | empirical hard cap 3.6880704671415807e-6 | finite, measured spatial diagnostic |
| pressure projection fraction | not separate | hard <=1e-6 if strong >=1e-12 |
| relative strong divergence | diagnostic | diagnostic |
| BC/gauge/material/mass/resolution/energy/work | hard | unchanged hard |

The old cap is retired **only for new A9 acceptance**, because pressure-space
projection and predeclared spatial convergence both qualified its diagnostic
interpretation. It was not increased by a percentage. Its value, source and
historical STOP remain stored with active_hard_gate=false. No exact pointwise
incompressibility is claimed. All four new moderate roots were solved fresh.

## 9. Two source freezes and provenance

Execution base HEAD: `5f2d9498e1d0715f57b125ee49b4ddfd533467ae`.
Audit source SHA: `6f84f21df52f5811e0095814b45c55df5b565a8a8bfe95f67489238c8d2b90a7`;
audit archive SHA: `c4e21e4908ba655ba70d99d8783267ec37a53d34e5b16d6418e8f11536983408`;
decision-rule SHA: `b2a4a12c7557647e437cb63df2689322136f6bbff5f66c3d48b265a87426a8ca`.

Only after audit PASS, new production runner files were added. Every archived
audit implementation file remains byte-identical. Production source SHA:
`64250f222cb8fe4eda98e2abd7ec0f91be3d0045bd65a45e81f098833795b717`;
archive SHA: `f81b460f8a73f7cddbcd0f77e6427601b207419abb39bfa48d55568138ad9f6e`;
new policy SHA: `557896045044c8af06efe09744944a0549160051e64cfde5e02b808cbae9935d`.
All six core files match STEP3A8. DOLFINx0.10.0 / PETSc3.24.0 / MUMPS5.8.1,
float64, one rank, identical pinned image. Final git commit is not execution HEAD.

## 10. Moderate full bridge

P1 atol=rtol1e-12; P2 atol=rtol3e-13. All other controls unchanged. This is
nonlinear-root equivalence, **not** temporal qualification at dt1e-6.

| metric | P1 abs | P1 rate | P2 abs | P2 rate |
| --- | --- | --- | --- | --- |
| SNES reason / iterations | 2 / 2 | 2 / 2 | 2 / 2 | 2 / 3 |
| initial residual | 0.784389086476 | 0.00242370482416 | 0.784389086476 | 0.00242370482416 |
| final residual | 1.08171994744e-13 | 9.00383837736e-13 | 1.08171994744e-13 | 9.69539175721e-15 |
| u L2 | 7.563407043299e-9 | 7.563407043201e-9 | 7.563407043299e-9 | 7.563407043299e-9 |
| pi L2 | .001091516415835943 | .001091516415840941 | .001091516415835943 | .001091516415835934 |
| phi L2 | .8067374409515263 | .8067374409515263 | .8067374409515263 | .8067374409515263 |
| mu L2 | .15267139920381542 | .1526713992038782 | .15267139920381542 | .15267139920381542 |
| strong divergence | 3.814612027896e-6 | 3.814612027851e-6 | 3.814612027896e-6 | 3.814612027896e-6 |
| projected divergence | 1.178343294066e-21 | 8.387206494785e-22 | 1.178343294066e-21 | 5.007041365622e-22 |
| weak continuity Riesz | 1.178343294066e-21 | 8.387206494785e-22 | 1.178343294066e-21 | 5.007041365625e-22 |
| f_Q | 3.08903e-16 | 2.19870e-16 | 3.08903e-16 | 1.31260e-16 |
| mass [m²] | -.6060961027499453 | -.6060961027499453 | -.6060961027499453 | -.6060961027499453 |
| D_CH [W/m] | .003669671786703607 | .003669671786710035 | .003669671786703607 | .003669671786703464 |
| energy [J/m] | .06564277766876435 | .06564277766876439 | .06564277766876435 | .06564277766876435 |
| certified cells | 8.328131075834218 | 8.328131075834218 | 8.328131075834218 | 8.328131075834218 |

Differences below are norms of differences, not differences of norms.

| cross-formulation metric | P1 | P2 | threshold |
| --- | --- | --- | --- |
| relative phi L2 | 5.119489756612193e-15 | 3.831096269576499e-17 | 1e-10 |
| relative mu L2 | 2.797916446095754e-12 | 2.441063827453512e-15 | 1e-10 |
| relative D_CH | 1.751487843648647e-12 | 3.907985046680551e-14 | 1e-10 |
| scaled u difference | 5.96841655906335e-17 | 1.371588869658539e-20 | 1e-9 |
| scaled pi difference | 2.817053985087031e-13 | 3.734686447750502e-16 | 1e-9 |
| energy absolute | 4.163336342344337e-17 | 0 measured | 1e-12 |
| mass absolute | 0 measured | 0 measured | 1e-12 |
| relative strong-divergence difference | 1.185444918327219e-11 | 1.654269270609195e-14 | 1e-4 |

Every comparison and each independent physical candidate passed. Absolute
P1/P2 strong divergence is identical; rate P1/P2 change is approximately
1.18e-11 relative (diagnostic). Initial physical guesses differ as designed.
Core regression remains residual <=1.27780e-15, Jacobian action 4.27298e-17,
final finite-difference action error 2.83560e-13; no row rescaling was introduced.

## 11. Exact target tiny step

dt=4.657657028534679e-12. Actual atol1e-10, rtol1e-9, stol0, max_it30,
newtonls/bt, preonly/LU/MUMPS. Accepted, reason2, one Newton iteration;
initial/final SNES norms 9.012498970230127e-4 / 1.0183543158505285e-14.
True retained rate is archived directly, separately from reconstructed phi.

| quantity | full rate result |
| --- | --- |
| phi L2 difference from isolated | 0 measured |
| mu L2 difference from isolated | 4.292793196705405e-18 |
| E_kin [J/m] | 4.32524656805269e-25 |
| D_CH [W/m] | .11187268676021649 |
| D_visc / D_slip [W/m] | 8.427334708945781e-19 / 8.846738171521091e-30 |
| max speed at DOFs [m/s] | 3.886238055847416e-11 |
| strong / orthogonal divergence | 6.099865735412659e-10 / 6.099865735412658e-10 |
| projected / weak divergence | 1.419221974328787e-25 / 1.419221974328801e-25 |
| projection fraction | 2.326644611355164e-16 |
| mass change/domain | 1.387778780781446e-15 |
| certified cells / wall crossings | 8.328131075834218 / 2 |
| total-energy subtraction Delta E [J/m] | -5.23719956291302e-13 |
| stable work-evaluator Delta E [J/m] | -5.235196577650138e-13 |
| B_BE [J/m], not physical heat | 2.455017107779116e-15 |
| weak-work defect [J/m] | -3.942429374703267e-20 |
| decomposition roundoff [J/m] | 4.556839936383097e-21 |
| advective/diffusive phase Riesz ratio | 1.602368667303663e-12 |

All material, finite, mass, BC/gauge, topology, resolution, continuity and work
checks passed. This is one-step evidence, not transfer qualification.

## 12. Pilot, restart and first dt boundary

The primary dataset is `full_coupling_L0/`, one continuous prefix from t=0.
The first 50 intervals were reused, not rerun. A new-process diagnostic fork
from checkpoint25 reproduced steps26–28 with zero measured differences in
u/pi/phi/mu/rate, all four divergence quantities, integrated CH/viscous/slip
dissipation, energy and budget. The fork is not a replacement trajectory.

| pilot metric | result |
| --- | --- |
| accepted primary intervals | 127 |
| final pilot time [s] | 6.148107277665776e-10 |
| first block end / first post-boundary | 122 / 123, read from schedule |
| boundary crossed plus five intervals | PASS |
| first post-boundary dt [s] | 9.315314057069358e-12 |
| Newton iterations before / after boundary | 1 / 1 |
| maximum Newton iterations over prefix | 1 |
| previous rate guess unscaled | verified SHA match |
| restart | PASS, all measured comparison errors zero |
| rejected primary intervals | 0 observed over 127; not a full-L0 claim |
| max mass/domain | 1.9737298215558338e-14 |
| min certified cells | 8.328131075834218 |
| max absolute weak-work defect [J/m] | 8.936723280342688e-20 |
| continuous prefix sessions wall [s] | 753.6905173360065 |
| mean accepted interval wall [s] | 5.785959464960374 |
| p95 accepted interval wall [s] | 12.904329958901508 |
| full-L0 projected wall [s] | 15942.5699042599 — FAIL cost gate |

All every-step physical gates pass. The stored trajectory execution status
`incomplete` refers to full 1068-step L0; the requested boundary pilot itself
completed successfully. Session files are closed with status=complete.

## 13. Cost policy

Audit scientific wall=111.35227182901144s including initialization/equilibrium
and rejected-state diagnostics, below 3600s. Preliminary cap1800s, full-L0
cap14400s and total19800s are unchanged. Pilot forecast uses measured full
accepted-step timing with 1.25*max(mean,p95). Pilot is counted once in total.

| wall-time quantity | measured/projected seconds |
| --- | --- |
| audit, including preparations and read-only candidate audit | 111.35227182901144 |
| moderate/tiny/pilot/restart preliminary sessions | 881.0210211020094 |
| of those, continuous L0 prefix sessions | 753.6905173360065 |
| total measured scientific sessions, no double counting | 992.3732929310208 |
| conservative accepted-step forecast rate | 16.130412448626885 |
| remaining 941 intervals plus initialization forecast | 15188.879386923894 |
| projected complete L0 | 15942.5699042599 |
| frozen complete-L0 cap | 14400 |
| projected total including audit/preliminaries | 16181.252679854915 |
| actual complete-L0 runtime | NOT RUN |

The full-L0-specific cap fails although the total allowance would pass. Mean
alone cannot replace the frozen p95/safety-factor policy. Mean components:
SNES3.47736s, basic diagnostics0.16779s, resolution0.09825s,
continuity/projection observer0.14748s, topology0.01284s, BE work0.13791s,
archival1.69296s per accepted interval. Sparse field/checkpoint work makes the
timing distribution nonuniform; it was included, not discarded from p95.
No source/cache/diagnostics optimization was introduced after seeing this gate.
General pytest and final read-only report rendering are outside scientific
session clocks; these figures are wall time, not CPU time.

## 14. Full-L0 physical and strong-divergence results

**INCOMPLETE: 127/1068**, final accepted time6.148107277665776e-10 rather than
1e-4. No full-horizon extrema or completion claim is inferred. The authorization
check raised CostStop before any interval128 solve; no failed PDE attempt is
being retried. Latest accepted checkpoint127 remains intact.

The following statistics include t=0 and the 127 accepted endpoints ONLY.
The measured zero minimum is the actual initial u=0, not missing data.

| prefix quantity | min | mean | p50 | p95 | max |
| --- | --- | --- | --- | --- | --- |
| strong divergence L2 | 0 | 3.159098992251e-8 | 3.263664111388e-8 | 5.615297465931e-8 | 6.044641155439e-8 |
| projected divergence L2 | 0 | 7.008374459993e-24 | 6.224171746296e-24 | 1.671141214594e-23 | 2.284066203706e-23 |
| orthogonal divergence L2 | 0 | 3.159098992251e-8 | 3.263664111388e-8 | 5.615297465931e-8 | 6.044641155439e-8 |
| projection fraction | 0 | 2.140404327438e-16 | 1.883141617997e-16 | 3.943468526752e-16 | 4.979369080447e-16 |

Max weak continuity Riesz=2.2840662037059463e-23. Max strong divergence is at
step127: t6.148107277665776e-10, dt9.315314057069358e-12,
velocity L2=1.0062411340047978e-10, E_kin=5.0626060988162314e-21 J/m,
D_visc=7.890139235800417e-15 W/m, D_CH=0.04663873343088977 W/m.
All full-horizon min/percentiles/max remain NOT RUN. Weak continuity is
satisfied to measured tolerance; strong divergence is explicitly nonzero.

## 15. Coupling and energy closure

NOT QUALIFIED DUE TO COST STOP. Reverified isolated L0–L2 references: max common phi L2
3.691188008762449e-8; integrated-D_CH gap 7.033312454099817e-11 J/m;
final F gap 6.0652038946784614e-12 J/m. Initial excess remains
5.542888467657825e-8 J/m; frozen energy evaluation floor1e-13 J/m.
The full coupling comparison must use identical L0 times, no interpolation.
Combined proxies add measured coupling and isolated temporal discrepancy;
they are practical conservative proxies, not rigorous continuum error bounds.

Read-only **prefix-only** comparisons use the same scalar clock and the
predeclared common field checkpoints that have already been reached. No
full-horizon transfer gate is declared passed from these values.

| observable | measured prefix coupling | isolated full-horizon temporal reference | full transfer pass? |
| --- | --- | --- | --- |
| max common phi L2 | 1.7771940626594783e-17 | 3.691188008762449e-8 | NOT QUALIFIED: incomplete |
| integrated D_CH, at prefix endpoint | 7.260771697140268e-22 J/m | 7.033312454099817e-11 J/m at 1e-4 | NOT RUN at common final horizon |
| F_CH, at prefix endpoint | 0 measured | 6.0652038946784614e-12 J/m at 1e-4 | NOT RUN at common final horizon |

The full-horizon references in this table are context, **not** a valid comparison
of final errors at different times. On the same prefix window the max isolated
phi0–phi2 gap is 1.514170751460509e-10; its integrated-D gap at prefix end is
4.7398431867882595e-14 J/m. Practical prefix-only combined proxies are
1.514170806623663e-10 (max common phi) and 4.7398432593959765e-14 J/m (integral).
The unknown full-horizon coupling/combined discrepancies are NOT RUN.

| prefix-only diagnostic | measured max/final | unchanged full-horizon threshold |
| --- | --- | --- |
| phi coupling / initial disturbance | 1.831034439021824e-13 | 1e-3 |
| F_CH coupling / initial excess | 2.503710455079495e-10 | 1e-3 |
| significant D_CH relative coupling | 5.376077850158562e-11 | 1e-2 |
| E_kin / initial excess | 9.133516087065452e-14 | 1e-4 |
| hydroD / D_CH | 1.691756755704856e-13 | 1e-4 |
| cumulative hydroD / initial excess | 3.2856697946299e-17 | 1e-4 |
| advective/diffusive phase Riesz ratio | 4.387866205343724e-10 | diagnostic only |

| prefix energy quantity | measured |
| --- | --- |
| endpoint CH free energy [J/m] | .06564278264512508 |
| Delta E_total [J/m] | -3.942120241351432e-11 |
| integral D_CH [J/m] | 3.9501616171368146e-11 |
| integral D_visc [J/m] | 1.8212101212920134e-24 |
| integral D_slip [J/m] | 2.6559264407930295e-35 |
| total energy-budget defect [J/m] | 8.041375785564643e-14 |
| defect / actual prefix energy change | .00203986060628326 (0.204%; diagnostic, not final closure qualification) |
| sum B_BE [J/m], not heat | 7.912223485457362e-14 |

No t=0 power sample or first interval is omitted. Integrated dissipation is
trapezoidal physical power, not substituted by minus energy change.
Prefix plots and metrics are under [analysis/](validation_results/step3a9/analysis/).
All full-CHNS curves there stop at the measured prefix; none extends to 1e-4.

## 16. Tests and scope

Audit preflight: pinned19 passed; production-policy preflight: pinned25 passed,
including new full-trajectory weak/nonfinite rejection, material/SNES rollback,
and new-process checkpoint/restart across a dt change. Final local host suite:
**293 passed, 46 skipped, 2 deselected**, 17.93s. Full pinned STEP3–STEP3A9
suite: **224 passed**, one expected under-resolution warning from an intentional
negative fixture, 77.33s. JUnit files are retained in `step3a9/tests/`.
No remote CI has been observed for this unpushed source. Production mesh studies
and L0 are excluded from GitHub Actions. CI includes the new tiny policy tests.

No full CHNS L1/L2, isolated rerun, different FE family, unequal-density AGG,
force, BDF2, moving contact line, contact-angle quench, sensitivity series,
falling film or tank benchmark is run or qualified here.

## 17. Scientific verdict and NOT RUN

TAYLOR-HOOD DISCRETE INCOMPRESSIBILITY POLICY QUALIFIED BY PRESSURE-SPACE
PROJECTION AND STRONG-DIVERGENCE MESH REFINEMENT;

FULL CHNS PHASE-RATE BRIDGE PASS;
COUPLING TRANSFER INCOMPLETE DUE TO COST GATE.

NOT RUN: intervals128–1068, complete full L0, final-time full-vs-isolated
coupling and transfer-dominance gates, full-horizon kinetic/hydro extrema,
full-horizon energy closure and full combined uncertainty proxies. No claim
that isolated temporal qualification transfers to full CHNS has been made.
The successful divergence audit is not substituted for that missing proof.

MODEL NOT YET VALIDATED.
