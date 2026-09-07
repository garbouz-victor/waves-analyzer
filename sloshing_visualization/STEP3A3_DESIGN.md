# STEP 3A.3 — policy fixed before new PDE runs

Base HEAD: `29d40b6cbddb882d5b6248dc2f52df5f8f93b2eb`.
Verdict remains **MODEL NOT YET VALIDATED** until the new measured gates pass.
All STEP 3/3A.1/3A.2 reports, data and prepared coefficients are read-only.

## Invariants and sequence

Use the STEP 3A.2 `prepared_symmetric` 96x48 crossed P2 equilibrium, its
degree-12 quadrature, and exactly the old delta=1e-3 compact-bump perturbation.
No physical coefficient, equation, wall condition, pressure/AGG term, FE degree,
equilibrium coefficient, amplitude or shape changes. No clipping, smoothing,
post-solve mass/energy correction, t=0 omission, or modified physical heat.

The old artifacts contain the perturbation and initial-phi SHA256, but not a
standalone saved psi coefficient array. Recover those bytes by replaying the
unchanged STEP 3A.2 deterministic initializer, require BOTH recorded hashes to
match exactly, and freeze the recovered arrays under step3a3 before analysis.
This is recovery of the same input, not design of a new perturbation. A hash
mismatch is a STOP. The old initializer source is in its saved source archive.

Sequence: input/operator checks -> relevant spectrum and converged Krylov
reference -> isolated CH low-mode benchmark -> linear-verified immutable plan
-> full nonlinear CHNS plan refinement. A failed prerequisite prevents all
downstream scientific runs. Tiny CI cases are separate from scientific runs.

## Operator and spectral policy

Assemble sparse M0, K and the UFL second variation H of the SAME free energy.
L = -mobility M0^-1 K M0^-1 H; never form dense production L or its inverse.
Mass projection Pq=q-1*(m^Tq)/(m^T1) is analysis-only. At least 20 deterministic
random vectors check conservation and bilinear symmetry.

* Relative mass defect: |m^T Lq|/(sqrt(area)*||Lq||_M) <= 1e-11.
* Relative ||PLPq-Lq||_M/||Lq||_M <= 1e-11 for projected q.
* M/K/H symmetry: Frobenius skew/size and bilinear defects <= 1e-11;
  bilinear normalization uses ||q||*||A r||+||r||*||A q||.
* Hessian central finite difference uses h=1e-3,5e-4,2.5e-4,1.25e-4;
  relative M-dual norm defect <= 1e-6 with decreasing errors until roundoff.
* Nonlinear semidiscrete rate central difference at the same amplitudes must
  approach L psi, final relative L2 error <= 1e-6. Also report the one-sided
  rate/delta sequence and mass defects; its O(delta) term is not mislabelled
  a central-difference failure. No time advance is used in these tests.
* Relative linear energy identity defect <= 1e-10 for the full operator.
* Early hydrodynamic power/CH power <= 1e-4 and E_kin/initial excess <= 1e-4
  justify a CH-only planner; otherwise stop for coupled-spectrum investigation.

Own M-inner-product Arnoldi uses two-pass modified Gram-Schmidt and mass
projection, starting at m=32/48/64. Increase to 96/128/192/256/384/512/768/1024
only if needed, without changing the 1e-6 self-convergence gate. Record all
levels, Ritz residuals and ||V^T M V-I||; orthogonality target <= 1e-11.
The spectral investigation resource guard is 1 hour and 2 GiB basis storage.
Failure at the finite guard is a documented blocker, never a reference PASS.

Ritz residual is ||Lv-lambda v||_M/(max(|lambda|,1 s^-1)*||v||_M).
A converged Ritz mode requires residual <= 1e-6 and relative eigenvalue change
<= 1e-6 between the final two levels. Seed significance is its reconstructed
modal L2 contribution / ||q0||_M; significant means >=1e-8, or its initial
quadratic dissipative self-contribution >=1e-8 of D2(q0). Report the complete
Ritz set, including unconverged/insignificant entries; do not discard failures.
For converged significant modes require Re(lambda) <= 1e-8*max(|lambda|,1 s^-1)
and |Im(lambda)|/max(|Re(lambda)|,1 s^-1) <= 1e-8. Significant unstable/complex
evidence triggers investigation/STOP, not a new timestep.

Reference horizon t_end=1e-4 s. Times include zero and a dense logarithmic
grid from 1e-3 of the fastest relevant tau through the entire horizon.
Last-two Krylov levels: max ||q_m-q_prev||_M/||q_m||_M <= 1e-6 (a 1e-14
initial-norm floor only avoids zero/underflow). Also require relative E2 and
D2 self-differences <= 1e-6 and record projected energy identity defects.
E2=q^T H q/2, r=M0^-1 Hq, D2=mobility*r^T K r. Compare E2/D2 with the old
nonlinear tau_E=4.862048987156724e-7 s; relative agreement target <=1e-2.

## Isolated CH low-mode policy

Select a real negative converged eigenmode nearest tau=3e-4 s, preferring
1e-4..1e-3 s, excluding a translation-like near-zero mode. If ordinary seeded
Arnoldi does not resolve such a mode, a sparse shift-invert eigensolve may
independently refine it; verify its FULL operator residual <=1e-8.
Normalize in M and set max |delta_mode*q_mode|=1e-4 before running. Additional
amplitude-halving checks use 5e-5 and 2.5e-5, without changing the stiff bump.
Horizon is one modal tau, BE dt=tau/20,/40,/80. Compare nonlinear isolated CH,
linear BE factor (1-lambda*dt)^-n and exp(lambda*t).
The two amplitude-halving checks use the same tau/20 schedule as the coarsest
low-mode run; their comparator is the exact linear BE law, not exp(lambda*t).

Require mass/domain <=1e-10, >=8 certified cells, energy growth <=1e-9 J/m,
correct decay, decreasing modal errors with observed order >=0.8, decreasing
integrated-D differences with last relative <=5%, final continuum budget/
actual |DeltaE| <=5%, and nonlinear-vs-linear BE modal error <=1e-3 at the
smallest amplitude with decreasing amplitude error until roundoff. A tiny
neutral-wall CI case checks the same algorithm, not production qualification.

## Immutable spectral plan

z_target=0.05. Initial candidate dt0=z_target/rho_relevant, where relevant
modes are converged and significantly seeded as above. Do not use the full
mesh spectral radius blindly. Block dt may increase by at most a factor two.
Vanished-content heuristic cutoff is 1e-4 of initial significant content;
it is only a proposal, never acceptance. Candidate schedules are checked by
linear BE against the converged exponential reference before nonlinear use.

Required predicted errors: endpoint and max checkpoint relative L2 <=1e-3;
max checkpoint quadratic energy relative error <=1e-3; integrated D2 relative
error <=1e-2. Physical power uses trapezoids including t=0. No dense inverse.
Candidate refinement/rejection is linear analysis only. Save all rejected
candidates; a final failed candidate forbids nonlinear execution.

Freeze every block boundary, dt and integer interval count in a deterministic
JSON plan with operator/equilibrium/input/Krylov/policy/plan SHA256. Save the
full input provenance, predicted errors, z values and single-step BE modal
errors. No observed nonlinear field enters schedule construction or execution.
If SNES fails, record failure, never alter dt.

Cost uses measured STEP3A.2 seconds per accepted step (including diagnostics).
Pre-run guard: predicted level-0 runtime <=3 hours and total requested series
<=10 hours, with available memory checked. An exceeded guard is STOP/report,
not permission to loosen any numerical tolerance or launch blindly.

## Full nonlinear CHNS qualification

Only after every preceding gate passes: t_end=1e-4 s, level0=plan,
level1=all block dt/2, identical boundaries; level2=all block dt/4 if the first
two are insufficient and the resource guard allows it. Three levels are needed
to measure decreasing successive field/integral differences, so run level2
when that evidence is not yet available. No alternative block layout.

All levels must have coefficient-identical initial psi and phi, mass/domain
<=1e-10, >=8 certified cells at EVERY accepted state, and no physical energy
growth beyond the existing 1e-9 J/m allowance. Require decreasing final budget
defects, decreasing successive D_CH and endpoint L2 differences, last relative
D_CH difference <=5%, final budget/actual |DeltaE| <=5%.

Record BE work decomposition on intervals 1..20 and logarithmically spaced
later indices: DeltaE, D_old/new, both time integrals, B_BE, weak-work defect,
continuum defect, decomposition roundoff. Early absolute weak-work target
<=1e-12 J/m and full algebraic work balance consistent with configured SNES
residual; decomposition roundoff <=1e-12 J/m. Local first-interval ratios are
diagnostic, not a new hard gate. Physical heat is D_CH+D_visc+D_slip ONLY;
B_BE is a signed numerical work remainder and is never added to heat.

Report early linear/nonlinear L2/H1, excess-energy and power comparisons without
requiring exact equality. Preserve every accepted state diagnostic and enough
coefficient snapshots for all common endpoint/reference comparisons.

Maximum possible success verdict:
**BASIC NONZERO CH TRANSIENT TEMPORALLY QUALIFIED; MOVING-CONTACT STEP 3A
BENCHMARKS STILL REQUIRED**. No new spreading, angle sweep, film, tank, gravity,
water-air or unequal-density experiment is authorized in this iteration.

## Linear-reference implementation investigation (before any nonlinear run)

The unrestarted polynomial Arnoldi reference is being tested through m=1024;
already at m=768 the full-horizon L2 difference is 0.240, despite machine-level
M orthogonality. The measured upper Ritz rates are approximately 1.07e10 s^-1,
so t_end*rho is approximately 1.07e6. A single polynomial space is inefficient
for the long end of this interval. These failed levels will not be overwritten.

The installed SLEPc 3.24 MFN Krylov exponential supports restarts and a BV
mass-matrix inner product. As an additional LINEAR investigation, test this
exponential action with M inner product, always-reorthogonalized modified
Gram-Schmidt, and ncv=32/48/64 (96/128 if necessary), internal tolerance 1e-11.
Require the SAME <=1e-6 L2/energy/power self-convergence, plus independent tiny
dense-exponential tests. Keep the one-hour resource guard for this bounded
alternative. It is not nonlinear BE, a modified model, a looser gate, or a
replacement of physical dissipation by energy change. Unrestarted Ritz results
remain separately labelled; no nonlinear stage may proceed on a failed reference.

The mandatory tiny nonidentity-M dense-exponential regression rejected the
direct MFN/BV matrix setting (relative M-norm error approximately 0.1775): MFN's
initial normalization is Euclidean. That adapter is NOT used for scientific
data. The independent restarted MFN alternative therefore retains its normal
Euclidean internal basis and is certified externally in FEM M/L2 norm and
the exact H/K energies/powers. The primary explicit M-Arnoldi Ritz analysis
is unchanged. This changes no operator, no scientific tolerance and no PDE.

MFN's restarted polynomial action also becomes expensive at the long end
(ncv=32: 876.77 s to t=4.5854e-5, growing projected work). Before any nonlinear
run, permit a second independent reference construction: mass-orthogonal
rational Krylov, with fixed positive resolvent times log-spaced from 1/rho_Ritz
to 1e-4 s, repeated cyclically. The resolvent is the same sparse coupled
linear BE solve used only as a BASIS GENERATOR, not as the reference evolution.
Evolution is still exp(t*V^T M L V) and must pass the unchanged 1e-6 L2/E2/D2
self-convergence checks at m=32/48/64 (96/128 if needed). Retain full-operator
Ritz residuals; verify the implementation against a tiny dense exponential.
Use a bounded one-hour analysis and at most 2 GiB of cached resolvent factors.
No nonlinear execution is authorized merely by this alternative existing.

The first unbalanced sparse resolvent prototype reached dimension 8 in 215.49 s
and 4.19 GiB RSS. It is stopped for a solver-conditioning/cost improvement:
use r_tilde=sqrt(dt)*r in the coupled linear system, producing blocks
[[M,sqrt(dt)*mobility*K],[-sqrt(dt)*H,M]]. Eliminating r_tilde gives the SAME
linear BE resolvent. The original operator, M norm, reference exponential,
initial coefficients and thresholds do not change. A dense tiny comparison
checks the equivalence before a new versioned rational-reference attempt.

The balanced resolvent uses the already pinned PETSc/MUMPS sparse factorization,
with at most three simultaneous cached factors (destroyed between batches).
Its tiny dense comparison includes dt=1e-10, 1e-3 and 1e-1. No PETSc/SLEPc
package installation or model change is involved.
