# STEP 3A.2 — exact P2 certification and prepared variational equilibrium

**MODEL NOT YET VALIDATED**

Discrete equilibrium preparation, exact P2 spatial certification, and both BE
and BDF2 equilibrium preservation **PASS**. The prescribed three-level nonzero
CH relaxation series **does not pass energy/dissipation qualification**. The
scientific sequence stops there. No new scientific contact-angle quench,
spreading, tank, film, water–air or gravity series was run.

Base HEAD: `04904c00c90b54da0ea03525399407ef2ce98b6c`. All new measurements are
under [validation_results/step3a2](validation_results/step3a2/summary.json).
Historical STEP 3 and STEP 3A.1 results are unchanged. The policy was recorded
before the corresponding runs in [STEP3A2_DESIGN.md](STEP3A2_DESIGN.md).

## 1. Goal

Establish a certified spatial diagnostic, a stationary discrete mass-constrained
state of the existing free energy, preservation by the existing transient
discretization, and then qualification of a controlled nonzero relaxation.
The weak form, AGG structure, pressure transform, force sign, sigma normalization,
bulk/wall energy, angle convention, slip, mobility and material parameters were
not changed. No clipping, solution mass correction, smoothing, energy correction,
modified dissipation, t=0 removal, wall kinetics or hysteresis was introduced.

## 2. Why geometric 60° was not variational equilibrium

The historical `compatible_60_degree_energy_baseline` is **geometrically
compatible, variationally unprepared**. Its zero-contour angle is 60°, but

`L = (3 sigma/4) [R cos(theta_initial)/r - cos(theta_e)] (1-phi²)`.

For equal initial/equilibrium angles, the bracket vanishes at r=R, not throughout
the diffuse wall trace. A geometric circle therefore does not solve the natural
wall variation or the complete discrete chemical equation. Historical results
and their names inside archived records have been retained, not rewritten.
Their measured 577.494578 W/m initial CH power remains an unqualified startup.

## 3. Exact P2 polynomial certification

[p2_certification.py](src/sloshing/multiphase/p2_certification.py) reconstructs
`a xi² + b xi eta + c eta² + d xi + e eta + f` from the actual six nodal
coefficients. Coefficients come directly from the FE dofmap in reference-node
order, avoiding basis-evaluation cancellation for constant fields.

Candidates are all three vertices, each interior edge critical point
`t=-B/(2A)`, and the interior solution of the 2×2 stationary system when its
Hessian is well conditioned. Singular Hessians need boundary extrema; stationary
lines reach the boundary. A nearly singular nonzero determinant invokes a
labelled `certified_conservative_fallback` Bernstein bound, never a regularized
solve. Finite-precision bounds include an explicitly recorded outward roundoff
allowance. Witness type, reference/physical location and unpadded value are saved.
CG1 retains exact vertex extrema.

Pure regressions cover constants, linear fields, convex/concave interior extrema,
edge extrema, saddles, near singularity, stationary lines, vertex/edge band
crossings with a pure centroid, and a genuine Bernstein false positive. Seeded
tests compare 350 random P2 polynomials with 5151 dense points each. Dense
evaluation is a test only, never the detection algorithm.

## 4. Certified normal-spacing method

Physical P2 gradients are affine and are reconstructed in extended precision.
Their three vertex values define a convex hull. If that hull contains or nearly
touches zero, the certificate uses the triangle diameter. Otherwise its complete
angular cone is found. For each of the six oriented vertex pairs, projected
width is maximized at a cone endpoint or at the direction parallel to that pair.
The largest candidate is the certified width, with outward numerical allowance.

The hard gate is `4.164065537… epsilon / h_cert >= 8`. The seven-point
directional estimate remains separate. Tests verify exact linear projection,
rotation covariance, origin-hull fallback, and 250 random hulls × 4000 sampled
directions without an underestimated width. No sampled value can supply the
hard certificate.

## 5. Historical resolution re-audit

The original t=0.5 checkpoint was loaded without advancing the PDE. Its SHA256
was identical before and after the audit.

| method | minimum transition cells | active cells |
|---|---:|---:|
| old centroid indicator | 7.4882791402 | historical indicator |
| STEP3A.1 Bernstein + sampled normal | 3.8659053283 | 3244 |
| exact P2 + certified complete cone | 3.865905328299158 | 3244 |

**Bernstein false positives: 0.** The problematic element is a true transition
cell, not a Bernstein artifact. Cell 8612 has actual range
`[-0.9496029932541654, -0.8575853314780866]`; its minimum is on an edge, maximum
at a vertex. `h_cert=0.026928139622531787 m`. There were no range or diameter
fallbacks among historical active cells. The old all-time verdict is unchanged.

Both the initial expression-evaluation audit and the final direct-coefficient
[audit](validation_results/step3a2/p2_resolution/historical_reaudit/actual_coefficients.json)
are retained. Their minima differ only in the last floating-point digits.

## 6. Discrete constrained free-energy problem

The functional is exactly

`F_h = integral [lambda/(4 epsilon)(phi²-1)² + lambda epsilon/2 |grad phi|²]
       + integral_wetting f_wall(phi)`.

The augmented equations are `dF_h(phi)[chi] - mu_star integral chi = 0` for
every phase test function and `integral phi = M_target`. The scalar is a single
global algebraic DOF. There is no cellwise P0 multiplier and no new wall law.
Quadrature degree is 12 for the stationary and transient forms and energies.

Parameters remain rho_l=rho_g=1 kg/m³, eta_l=eta_g=1 Pa s, sigma=.1 N/m,
epsilon=.025 m, M=1 m²/(Pa s), L_s=.01 m, theta_e=60°, g=a_x=0. Domain:
[-.6,.6] × [0,.6] m. This is the matched-density demonstration, not water–air.

## 7. Equilibrium solver algorithm

[equilibrium.py](src/sloshing/multiphase/equilibrium.py) uses serial PETSc SNES
on N_phi+1 unknowns, with the exact bordered Jacobian and MUMPS. Every iteration
records mu_star, mass, mass defect, SNES norm, Newton update norm and energy.
The optional bounded scalar-root utility is independently tested; the measured
solutions use the augmented solve, so `outer_root_history=[]` is explicit.

The first adaptive-mesh solve stalled and was rejected. An instrumented repeat
measured normalized Riesz defect 3.12218e-7 and a Hessian eigenvalue 6.5590e-12.
The corresponding mode overlaps horizontal translation by 0.9999923 in L2.
Changing Newton globalization to PETSc's [affine-covariant nleqerr line search](https://petsc.org/release/manualpages/SNES/SNESLINESEARCHNLEQERR/)
alone did not cure this adaptive-mesh stall. All three failed records remain.

Uniform crossed triangulations then solved the **full** unconstrained-in-position
equations in 8 and 9 iterations. No center pin, symmetry projection, phase shift
or mass adjustment was applied. Only the mesh was changed. The new mesh option
defaults to the original right-diagonal triangulation for historical configs.

## 8. Target mass

The radius .3 m analytic 60° tanh is interpolated in each actual P2 space and
integrated with the same degree-12 rule. Targets are:

| mesh rectangles | M_target [m²] | final mass minus target [m²] |
|---|---:|---:|
| 96×48 crossed | -0.6060961027495515 | -3.9490633e-13 |
| 108×54 crossed | -0.6060961027844581 | +2.4213964e-13 |

The first target corresponds to liquid inventory approximately .0569519486252 m².
It is not an analytic-area constraint. The tiny difference between meshes is
their discrete interpolation/integration difference. No post-solve correction
was made; normalized final mass errors are below 5.5e-13, versus the 1e-10 gate.

## 9. Prepared equilibrium result

The first fully certified state is
[prepared_symmetric](validation_results/step3a2/equilibrium60/prepared_symmetric/prepared.json).
It uses 18432 triangles, 37153 phase DOFs and 157973 transient unknowns.
Preparation took 11.03 s; the finer preparation took 14.35 s.

| metric | analytic 60° tanh | prepared equilibrium |
|---|---:|---:|
| phase mass [m²] | -0.6060961027495515 | -0.6060961027499464 |
| apparent angle, 2–10 epsilon [deg] | 59.99873902 | 60.02592815 |
| wall residual L2 [N/sqrt(m)] | 8.47548331e-4 | 6.86828041e-5 |
| wall residual max, sampled [N/m] | 2.19824085e-3 | 5.99325041e-4 |
| stationary Riesz L2 [Pa m] | 5.33355499e-2 | 4.73274932e-14 |
| normalized stationary Riesz | 1.88569645e-1 | 1.67327957e-13 |
| mu coefficient spread [Pa] | 6.24125754 | 0 |
| D_CH at t=0 [W/m] | 1024.046613755 | 4.23566083e-28 |
| free energy [J/m] | .0668511936194712 | .0656427272556616 |
| certified transition cells | 8.3281310758 | 8.3281310758 |

This comparison uses the **same new mesh**, so its analytic initial power is
1024.05 W/m, not the historical different-mesh value 577.49 W/m. The prepared
energy is lower by .00120846636381 J/m. `mu_star=.17992028183404365 Pa`.

## 10. Weak stationarity residual

The L2 mass-matrix solve defines `(r_h,chi)=EL_defect(chi)`. The scientific norm
is `||r_h||_L2 / [(sigma/L_reference) sqrt(area)]`, with L_reference=.3 m and
mu_scale=1/3 Pa. The threshold 1e-8 is unchanged. Both meshes achieve less than
2e-13. Raw SNES norms are recorded as solver diagnostics, not scientific norms.

Independent five-point differences of the **scalar energy**, in deterministic
zero-mass bulk/interface/wall directions, give normalized first variations
2.098e-11, 2.111e-11 and 1.881e-12. No Newton residual vector is reused for these
checks. All three symmetric energy second differences are positive; this is
directional local-minimum evidence, not a proof of a global minimum.

## 11. Wall residual

The strong-style trace uses actual P2 polynomial values and affine gradients.
Degree-12 edge quadrature integrates its squared norm; the maximum is labelled
as a 65-points-per-edge diagnostic. Natural BCs are not required to vanish
pointwise at machine precision.

Prepared wall L2 improves by a factor 12.34 over the analytic FE field. It drops
from 6.86828e-5 to 5.45025e-5 on refinement; sampled maxima drop from 5.99325e-4
to 5.08206e-4. This supports the variational preparation and its spatial behavior.

## 12. Apparent contact angle

The existing circle-fit method gives 60.025928°, 60.025928° and 60.037729° in
the 2–6, 2–10 and 3–8 epsilon windows. The first two windows select the same
eligible points because the cap height is below 6 epsilon. On the finer mesh,
the 2–10 result is 60.021237°. These are apparent diffuse-interface angles.

First-mesh contacts are (-.2421846498, +.2421846635) m; finer contacts are
(-.2422101075, +.2422101061) m. Finite-epsilon equilibrium is not forced to be
the original analytic radius-.3 cap. Variational stationarity is the criterion.

## 13. Equilibrium spatial resolution

| mesh | free energy [J/m] | mu_star [Pa] | certified cells | wall L2 |
|---|---:|---:|---:|---:|
| 96×48 crossed | .0656427272556616 | .1799202818340 | 8.3281310758 | 6.86828041e-5 |
| 108×54 crossed | .0656425800455409 | .1799193937265 | 9.3691474603 | 5.45025137e-5 |

There are no diameter or polynomial-range fallbacks among active cells of these
prepared states or their measured transient histories. Diameter fallback is
exercised by the singular/near-zero synthetic cases and under-resolved tiny
FEM tests. Two levels establish the reported differences, not full continuum
convergence. The first fully certified mesh is used for every transient level.

## 14. BE equilibrium-preservation test

`initialize_prepared_equilibrium` sets u=0, phi=phi_eq, mu=mu_star coefficientwise,
pi=0 and both history vectors equal to the current state. It does not reproject
mu. All 100 BE steps, including t=0, retain the state with zero measured drift.

| scheme | dt [s] | t_end [s] | max speed [m/s] | phi L2 drift | mass error/domain | max D_CH [W/m] | energy drift [J/m] | abs budget [J/m] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BE | 1e-5 | 1e-3 | 0 | 0 | 5.48481e-13 | 4.23566e-28 | 0 | 4.23566e-31 |
| BDF2 | 1e-5 | 1e-3 | 0 | 0 | 5.48481e-13 | 4.23566e-28 | 0 | 4.23566e-31 |

H1 drift, kinetic energy, viscous/slip power, contact displacement and all angle
drifts are also zero. These are measured zeros; each solve accepts its already
converged state with zero Newton iterations. Absolute and positive-energy-scale
budget gates pass. The near-zero physical change is not used as a preservation
denominator. Every accepted row, nonlinear log, resolution certificate and
geometric observation is saved. Runtime: 132.28 s for BE.

## 15. BDF2 preservation test

One ordinary BE interval supplies the history, then 99 constant-step BDF2
intervals follow. No microsteps are used. The same gates pass in 141.21 s.

A preliminary tiny test with unnecessarily tight transient atol=1e-13 hit a
BDF2 cancellation floor: its phase-block residual was 1.3991e-12 for identical
histories. Chemical blocks matched to sign, and momentum was negligible. The
scale analysis is recorded in `conditioning/tiny_bdf2_roundoff.json`. Scientific
runs use the original transient atol=1e-10/rtol=1e-9; scientific thresholds and
stationary atol=1e-13 were not relaxed. No time-discretization expression changed.

## 16. Smooth perturbation definition

After both preservation passes, `phi0=phi_eq+1e-3 psi_h`. Two C-infinity compact
bumps have centers (-.12,.32) and (.12,.43) m and radii (.09,.06) m. Inside
elliptic r²<1 each bump is `exp(1-1/(1-r²))`, and it is zero outside. Their
opposite weights use their discrete FE integrals; coefficients are normalized
before subtracting the residual discrete mean.

Integral psi=-7.13e-18 m²; actual phase-mass change=-5.06e-21 m². The subtracted
mean is -2.593e-17, so the wall trace is only roundoff-sized. Delta was never
changed. Exact initial extrema are [-.9790813279, 1.0156954720]; all matched
material coefficients are admissible. No clipping was applied. All three runs
have identical saved psi and initial-phi coefficient hashes. Nonconstant mu is
initialized with the existing consistent weak projection.

## 17. Temporal refinement of perturbed equilibrium

All levels use pure BE, the same prepared fingerprint and mesh, and t_end=1e-4 s.
Power integration includes every accepted step and t=0.

| dt [s] | integrated D_CH [J/m] | Delta E [J/m] | final budget defect [J/m] | max local defect [J/m] | end phi L2 difference to finest |
|---|---:|---:|---:|---:|---:|
| 1e-5 | 6.04869128e-7 | -4.29645692e-8 | 5.61904559e-7 | 5.58570663e-7 | 1.10912457e-6 |
| 5e-6 | 3.22656896e-7 | -4.31160014e-8 | 2.79540895e-7 | 2.77132746e-7 | 3.73651185e-7 |
| 2.5e-6 | 1.82024836e-7 | -4.31909640e-8 | 1.38833872e-7 | 1.37148139e-7 | 0 |

Successive endpoint L2 differences are 7.35620e-7 and 3.73651e-7. End energy
differences also decrease. Mass/domain remains below 6.11e-13, energy never
increases, and certified transition resolution is at least 8.32813 throughout.
Runtimes were 93.19, 121.29 and 243.74 s.

## 18. Energy/dissipation convergence

**FAIL.** Absolute budget defects decrease with observed orders 1.0073 and
1.0097, but the CH integral's last relative difference is **77.26%**, exceeding
the predeclared 5% limit. Final defects divided by actual energy decrease are
13.0783, 6.48346 and **3.21442**, versus a required <=.05. Error reduction is
measured; adequate integral convergence and closure are not established.

The initial power is .1140031391 W/m. Initial excess energy is 5.54289e-8 J/m;
their ratio is 4.86205e-7 s, below the smallest tested dt. At the coarsest level,
95.23% of the CH integral is its first trapezoid. This supports unresolved early
power integration as the dominant failure; the integral is not certified heat.

The entire physical energy change is below the historical significance floor
6.56428e-8 J/m. Consequently the unchanged legacy single-history relative-to-change
metric is N/A. The first run's original summary recorded that limited gate as
passed. It is preserved, and the separate `qualification_review.json` explicitly
rejects it under the predeclared **nonzero final actual-change** gate. That gate
was added to the runner before levels two and three; no scientific threshold
changed. The global summary cannot infer a pass from the historical floor.

The symmetric ratio `abs(Delta E+I_D)/(abs(Delta E)+I_D+1e-10)` is separately
saved with `diagnostic_only=true`. It does not replace an existing gate.

Two independent diagnostics sharpen the interpretation:

1. At perturbed t=0, solve `(phi_dot,s)=-M(grad mu,grad s)` without advancing.
   Differentiate the scalar free energy in that direction. The resulting
   `dF/dt+D_CH=-7.60503e-15 W/m`, relative defect 6.67089e-14. Initial chemical
   and CH variations are consistent with the free energy to roundoff.
2. From saved final BE intervals, the exact work identity is
   `Delta E + dt D_new + B_BE = weak_work_defect`. Here B_BE contains the kinetic
   increment norm, gradient increment norm, and exact quartic/cubic backward
   Taylor remainders. Its sign is not artificially constrained. The continuum
   trapezoidal defect decomposes as
   `weak_work_defect - B_BE + dt/2*(D_old-D_new)`.

For the finest **last** interval, B_BE=1.94211e-12 J/m, the power-quadrature gap
is 3.94263e-12 J/m, and weak-work defect is 8.27125e-20 J/m. Thus the local
continuum defect is about 2.00052e-12 J/m. This diagnostic separates BE work
and power quadrature at the end; it does not certify the unresolved initial
interval or replace the failed global continuum budget. Polynomial-remainder
identities have independent pure tests, including negative-remainder cases.

## 19. Remaining limitations

The required STOP is honored after the prescribed perturbation series. No
further smaller-dt PDE experiments or return to an analytic-circle baseline were
made. Read-only rate/work audits do not change any saved solution or energy.
No complete continuum, epsilon, mobility, slip, moving-contact, falling-film,
unequal-density or gravity validation is claimed.

The prepared checkpoint has its own schema
`step3a2-prepared-variational-equilibrium-v1`, phase coefficients, constant mu,
target mass, energy, stationarity diagnostics and integrity fingerprint. The two
successful H5 files are included with the results. Failed candidates are marked
incomplete and cannot initialize the transient solver. Reuse checks physical/FE/
quadrature/wall config, target mass, mesh/partition, algorithm source and content
integrity. Time scheme and dt are explicitly separate from equilibrium identity.
The augmented solver and prepared checkpoint currently require one MPI rank.

The successful equilibrium preparations and transient runs contain base git
SHA, source SHA256 and source archive, Docker image ID, configuration, mesh
fingerprint and solver history. Early failed preparations retain source
archives and Newton logs, but their failure summaries lack complete runtime
provenance; this metadata limitation is not silently filled retrospectively.
The source
hash identifies the working implementation beyond the base git SHA. Transient
metadata also records equilibrium_source_fingerprint. The first prepared
fingerprint is `df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc`.
Pinned stack: DOLFINx 0.10.0, PETSc 3.24.0, degree-12 quadrature, Docker image
`sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`.

Transient rank checkpoints, per-cell NPZ arrays and JIT caches remain local and
ignored by Git. All accepted diagnostics, geometric observations, nonlinear logs,
source archives, prepared coefficients, summary tables and PDFs are retained.
Ordinary pure CI and the STEP 3 Docker workflow include the new tests. Final
local test counts are recorded in
[TEST_RESULTS.md](validation_results/step3a2/TEST_RESULTS.md); remote GitHub Actions
was not executed because no push was performed. Existing tiny CI smoke cases
are distinct from new scientific spreading experiments.

## 20. Scientific verdict

**MODEL NOT YET VALIDATED**

Exact P2 detection and certified normal widths: PASS.
Mass-constrained equilibrium, independent first variations, spatial certification
and improved/refining wall residual: PASS. BE and BDF2 preservation: PASS.
Nonzero perturbation endpoint convergence and mass/resolution: PASS.
Integrated dissipation convergence and nonzero transient energy closure: FAIL.

The equilibrium incompatibility of the old initial state has been removed and
tested. Basic nonzero CH transient energy accounting is still unqualified.
STEP 3A.2 therefore does not earn its overall PASS, and moving-contact STEP 3A
tests remain required before any wider application.
