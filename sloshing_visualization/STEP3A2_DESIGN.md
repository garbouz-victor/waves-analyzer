# STEP 3A.2 policy recorded before measurements

Verdict remains MODEL NOT YET VALIDATED until all gates pass. Historical STEP 3
and STEP 3A.1 files are read-only. Their geometric 60-degree tanh is
geometrically compatible, variationally unprepared.

Physical equations, energy, wall law, material parameters and all historical
thresholds are unchanged. The hard spatial gate is the nominal 4.164065537
epsilon transition width divided by a certified upper bound on triangle width
over every possible affine P2 gradient direction: at least 8 cells. CG1 ranges
are vertex-exact. P2 ranges use vertices, edge stationary points and an interior
stationary point. Near-singular uncertainty uses labelled Bernstein bounds.
The width certificate uses the full gradient cone; an origin/near-origin hull
uses diameter. Sampled directions remain diagnostic only.

The prepared solve uses the same P2 space and degree-12 quadrature. Target mass
is the integral of the interpolated radius .3 m, theta=60 degree tanh on EACH
mesh. There is no post-solve mass adjustment. A true single algebraic scalar
mu is coupled to the phase coefficients and one global mass equation.

Predeclared equilibrium gates: abs(M-Mtarget)/area <= 1e-10; L2 Riesz
stationary residual divided by (sigma/L_reference)*sqrt(area) <= 1e-8;
constant mu coefficient spread <= 1e-12 Pa; D_CH(0) <= 1e-16 W/m; wall L2
residual at most half the interpolated analytic state's residual on the same
mesh. The wall maximum is a quadrature-sampled diagnostic, not a supremum
certificate. Mesh refinement must reduce the wall L2 diagnostic. Deterministic
zero-mass first variations must be <= 1e-8 after normalization by
(sigma/L_reference)*sqrt(area)*norm(psi)_L2. Energy finite differences are an
independent diagnostic (not the assembled Newton residual).

Preservation uses BE, then one ordinary BE history interval and constant-step
BDF2, initially dt=1e-5 s, t_end=1e-3 s (100 accepted steps). Gates:
mass/domain <= 1e-10; max speed <= 1e-8 m/s; relative phase L2 drift <= 1e-8;
energy drift / positive initial component sum <= 1e-8; absolute budget defect
<= 1e-10 J/m and defect / positive initial scale <= 1e-8; max D_CH <= 1e-16
W/m. Both t=0 and every accepted step are retained. H1, contact and angle
drifts are reported. No microstep startup. A failed preservation stops the
pipeline and triggers term-by-term residual inspection, not a smaller dt.

Only after preservation: phi0=phi_eq+1e-3 psi_h, two deterministic smooth
compact bulk bumps away from the bottom wall, normalized to max coefficient
amplitude 1. Their discrete integrals set opposite weights; remaining discrete
mean is subtracted before the transient. No solution correction. Material and
exact phase extrema are checked before the series. Consistent weak mu is used.
Pure BE refinement starts at dt=1e-5, 5e-6, 2.5e-6 s to the same t_end=1e-4 s.
Require decreasing final absolute budget defect, convergence of integrated CH
dissipation (successive difference decreases; last relative difference <= .05),
end phase difference decreases, mass/domain <= 1e-10, no energy growth above
the existing 1e-9 J/m allowance, and >=8 certified cells at every step.
Final budget/positive scale <= .01 and final budget/actual energy change <= .05
retain the established closure targets; the all-time relative-to-change metric
is also saved and assessed. No continuum gate is replaced by a BE endpoint sum.

The symmetric ratio abs(delta E+I_D)/(abs(delta E)+I_D+1e-10 J/m) is
diagnostic_only=true. Exact-equilibrium preservation uses absolute and positive
scale budget defects since its physical change is approximately zero.

Failure to obtain equilibrium: PREPARED EQUILIBRIUM NOT OBTAINED. Failure of
preservation or perturbation: stop the scientific sequence. Maximum permitted
success verdict: DISCRETE EQUILIBRIUM AND BASIC CH TRANSIENT VALIDATED;
MOVING-CONTACT STEP 3A TESTS STILL REQUIRED. No spreading, tank, film or
unequal-density/gravity runs belong to this iteration.

## Transient solver conditioning audit, before the scientific 60-degree runs

An initial tiny development test used transient SNES atol=1e-13. Term-by-term
inspection found a BDF2 phase-block algebraic norm 1.3991e-12 for identical
prepared histories at dt=1e-5; BE gives 2.04e-31. Chemical blocks agree with the
stationary form exactly (opposite sign), and momentum is 3.23e-32. Coefficients
are [150000, -199999.99999999997, 49999.99999999999]. Their weighted sum is
subject to cancellation of O(1/dt) terms. An 8-eps forward-error allowance times
sum(abs(a_i)) is 7.11e-10 per unit mass weight, so 1e-13 raw-vector stopping is
below the floating-point evaluation floor on this coarse mesh. This is not a
physical residual or evidence of field evolution.

The actual preservation and perturbation runs therefore retain the EXISTING
transient solver defaults atol=1e-10, rtol=1e-9. This restores the original
solver controls; all scientific mass, drift, energy, power and Riesz gates above
remain unchanged. Stationary augmented SNES retains atol=1e-13, rtol=1e-12.
The measured development audit is in conditioning/tiny_bdf2_roundoff.json.
No time step or weak-form expression is changed to address this roundoff.

## Stationary nonlinear conditioning, before the revised direct solve

The first residual-norm-backtracking solve and its instrumented reproduction
both stopped after 80 iterations. The normalized scientific Riesz defect is
3.12218e-7 (FAIL), not a prepared equilibrium. The phase Hessian has a near-zero
eigenvalue 6.5590e-12, and its eigenvector's L2 overlap with horizontal translation
is -0.9999923. This identifies a nearly neutral translation direction; it is not
resolved by loosening the residual gate or adding a center/mass correction.

The next direct solve uses PETSc's error-oriented affine-covariant `nleqerr`
line search for the same augmented equations and same Jacobian. This changes
only Newton globalization; all scientific and stationary SNES tolerances are
unchanged. A finite 80-iteration bound remains. The failed attempts and their
source snapshots remain in equilibrium60/prepared and conditioning_audit.
Algorithm reference: https://petsc.org/release/manualpages/SNES/SNESLINESEARCHNLEQERR/

The bounded nleqerr attempt also failed. An adaptive crossed mesh was audited:
only 4638/6010 triangles have reflected counterparts after local refinement.
The next spatial choice is uniform crossed triangles (no adaptive splitting),
96 by 48 rectangles, same domain and P2/quadrature/physical coefficients.
Its diameter .0125 m provides >=8.32813 nominal transition cells everywhere.
Cost estimate: about 158k transient unknowns, 342 MB CSR before LU fill;
the explicit guard is raised to 210k for this and a 108 by 54 equilibrium-only
mesh comparison. Host available memory was 8.4 GiB at review. Preservation
cost will be measured on its first accepted steps before the nonzero series.
No translation pinning, symmetry correction of phi, or new physical constraint
is applied: the full mass-constrained equations still determine every phi DOF.

P2 certification now reads reference nodal coefficients directly from the FE
dofmap, then computes physical affine gradients in extended precision. This
avoids spurious gradients from cancellation of basis-derivative sums for a
constant phase. Historical expression-based measurements are retained as a
versioned audit and compared with the coefficient-based certificate.
