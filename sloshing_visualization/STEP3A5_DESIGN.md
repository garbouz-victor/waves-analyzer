# STEP 3A.5 — policy frozen before new PDE runs

Base HEAD: `6106f8977be23a12446c2ea2d5aed0ea3311cf17`.
Overall verdict: **MODEL NOT YET VALIDATED**. STEP3–STEP3A4 reports and
scientific artifacts are historical/read-only. Only STEP3A5 receives new data.

## Scope and immutable inputs

Change only the algebraic BE phase unknown: r_phi replaces phi_new, with
phi_new = phi_old + dt*r_phi. Phase residual is (r_phi,s)+Mmob(grad mu,grad s).
Chemical residual is (mu,chi)-F'_h(phi_old+dt*r_phi)[chi]. The existing free
energy, exact wall variation, matched-density CHNS, sigma, epsilon, mobility,
slip, pressure transform, mesh, P2 spaces and quadrature 12 remain unchanged.
Old absolute solver and BDF2 are retained unchanged. No rate-form BDF2.

Require exact historical fingerprints:

- initial phi: `6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681`;
- psi: `831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0`;
- prepared equilibrium: `df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc`;
- frozen execution plan: `5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8`.

No clipping, smoothing, mass/energy correction, time-zero/interval omission,
schedule modification or formulation switch. No level0/half/quarter series
or coupling probe is authorized in this iteration. No contact-line, film,
tank, gravity, water-air or unequal-density experiments.

## Fixed nonlinear options and initial guesses

All primary runs: atol=1e-10, rtol=1e-9, stol=0, max_it=30;
SNES newtonls, backtracking (bt); direct preonly/LU/MUMPS, no variable/block
scaling or nonlinear preconditioning. Require converged reason >0, never
"almost converged". Log initial/final residual and both atol and rtol*r0.

First isolated rate guess: r0=-Mmob M0^-1 K mu_old, without mass projection;
mu guess is the unchanged consistent weak mu_old. Subsequent pilot guesses
use the previous ACCEPTED rate and mu, including block changes (no rescaling).
Zero-rate guess is permitted only in labelled pure/tiny tests, not a runtime
fallback. Record L2/max rate, mu norm and dt*max|rate| BEFORE first solve.

Reconstruct coefficientwise only after convergence. Keep rate unknowns separate
from physical checkpoints. Failed SNES never publishes phi, mu, time, histories,
accepted-step count or energy integrals. No retry with altered settings.

## Cheap equivalence and moderate-dt gates

Pure algebra: scalar exact/high-precision root and tiny-dt loss-of-digits tests.
Tiny FEM: assembled phase/chemical residual equality at representable moderate
dt <=1e-12 relative; J_rate=J_absolute*diag(dt,1), columns matched by scalar
DOF maps, <=1e-12 relative. Finite-difference Jacobian relative defect <=1e-6
at an appropriate central-difference step. No dense production Jacobian.

Production moderate comparison: ONE absolute and ONE rate step, dt=1e-6 s,
same frozen prepared/bump initial state and same SNES options. Require both
converged; relative phi L2 <=1e-10, relative mu L2 <=1e-10, relative D_CH
difference <=1e-10, absolute free-energy difference <=1e-12 J/m, absolute
phase-mass difference <=1e-12. Normalize phi/mu by the absolute solution norms;
D_CH by its absolute-run value. Failure stops the tiny primary experiment.

## Qualified sparse linear tiny step

Use historical qualified full sparse M0/K/H; BE dt exactly
4.657657028534679e-12 s. No dense full L/inverse. Verify linear phase-rate
system against existing balanced full sparse linear BE; endpoint relative L2
<=1e-10 and increment relative L2 <=1e-6 (increment denominator is much smaller).
Mass/domain <=1e-11. Keep seed q0=psi and physical amplitude delta=1e-3.

## Production isolated tiny-step proof

ONE step at dt=4.657657028534679e-12 s. Reproduce D_CH(0)=
0.11400313905309394 W/m to relative 1e-10; certified cells >=8 (historical
8.32813107583), two wall crossings, and exact input hashes before solving.
Require unchanged options, SNES reason >0, finite reconstructed state,
|Mnew-Mold|/area <=1e-10, certified cells >=8, exactly two wall crossings,
physical energy growth <=1e-9 J/m, |weak BE work defect| and decomposition
roundoff <=1e-12 J/m. Physical energy/work uses reconstructed absolute states.

Require relative L2 difference between nonlinear physical update/delta and
sparse linear update <=1e-3 (=delta), an explicitly O(delta) allowance.
Historical one-sided semidiscrete linearization error at delta=1e-3 was
9.8271e-6, so this predeclared gate allows nonlinear corrections without
claiming equality. Report both reconstructed-coefficient and dt*rate updates.

STOP after solving to reassemble/analyze, before authorizing pilot. On the
same accepted physical state record independent literal absolute UFL,
coefficient-difference, sparse mass-matrix and retained phase-rate residuals.
Report phase/chemical algebraic norms and L2 Riesz norms separately. Riesz
norms and one-coefficient-ULP/dt are diagnostics, not relabelled SNES norms.

Root-cause evidence additionally requires the historical unchanged absolute
form failed at this dt, the new form passes, physical/linear checks pass,
and literal absolute phase residual exceeds BOTH its rate counterpart by
at least 10x and the new SNES target. No threshold editing after measurement.
Compare historical REJECTED iterate only as diagnostic, never force a match.

## Pilot and unchanged cost gate

Only after all first-step/root-cause gates PASS: fresh 50-step isolated run
from identical t=0, first 50 exact level0 plan intervals, phase-rate throughout.
Every accepted state: SNES reason >0, finite/material checks, mass, energy,
physical powers, certified resolution and exact wall topology. BE work uses
STEP3A4 deterministic selections: 1..20 plus prescribed logarithmic/block-end
indices. Geometry at t=0, end and existing selected checkpoints. No adaptation.

Measure accepted-step and SNES costs, not the old failed attempt. Apply unchanged
STEP3A4 limits: isolated level0 <=2 h, pilot+three levels <=8 h; estimate using
max(pilot mean SNES, historical low-mode mean SNES), measured mandatory
diagnostics and sparse audits/geometry. Even a cost PASS does NOT authorize
the three-level series here. This first 50-step prefix does not cross a dt
block boundary; disclose that limitation, test dt changes only on tiny toys.

New bounded analysis stages <=1 h each and <=2 h total; no need to approach
these limits for one-step diagnostics. Full-CHNS one-step study <=2 h total.
Failure at any gate stops downstream PDE execution; no alternate solve path.

## Matched-density full-CHNS validation-only alternate path

Only after isolated proof/pilot review: add (u,pi,r_phi,mu) BE alternate solver
with phi_new=phi_old+dt*r_phi in advection, chemical, capillary and material
terms. Momentum time derivative and other original terms retain their BE
form; no velocity-rate reparameterization. Require matched density, g=a_x=0,
P2, quad12; do not generalize to unmatched density or BDF2. Publish physical
(u,pi,phi,mu) only after convergence, with original boundary/gauge conditions.

ONE moderate full-CHNS comparison at dt=1e-6: same phi/mu/energy/mass/D_CH
limits as isolated comparison. Velocity L2 difference <=1e-12 +1e-8*||u_abs||;
pressure L2 difference <=1e-12 +1e-8*||pi_abs|| (same gauge). Viscous/slip
power differences <=1e-18 W/m +1e-8 times absolute-run power. These solver-
scaled bounds avoid dividing by nearly zero hydrodynamic quantities.

After moderate PASS, ONE production full-CHNS rate step at exact tiny dt.
Same SNES/mass/resolution/finite/energy/BE-work gates. Compare isolated/full
physical fields and powers diagnostically; this single step is NOT a
coupling-negligibility qualification. STOP after that step and report.

## Provenance, monitoring, report and verdict

Store actual HEAD, source hash/archive, old/new formulation versions, mesh,
equilibrium/psi/initial phi/plan fingerprints, dt, options, PETSc/MUMPS stack.
SNES iteration records include residuals, available KSP reason/iterations,
Jacobian timing and line-search diagnostic trace; unavailable API quantities
are null with reason, not invented. Timing events are inclusive and separate
from application wall/CPU timers. No claimed exact huge condition number.

Maximum: TINY-dt BE ALGEBRAIC CONDITIONING BLOCKER RESOLVED BY PHASE-RATE
REPARAMETERIZATION; STEP3A NONLINEAR TEMPORAL QUALIFICATION STILL REQUIRED.
If pilot passes, additionally PHASE-RATE ISOLATED CH PILOT QUALIFIED;
FULL THREE-LEVEL SERIES STILL REQUIRED. Overall MODEL NOT YET VALIDATED.
If rate fails: RATE FORM DID NOT REMOVE TINY-dt SOLVER BLOCKER; no pilot.
