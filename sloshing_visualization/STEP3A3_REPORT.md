# STEP 3A.3 — CH spectral validation checkpoint

**MODEL NOT YET VALIDATED.** Operator, rational exponential reference and
isolated low-mode checks PASS. Full nonlinear stiff-bump qualification is
**NOT RUN**: the predeclared production cost gate stopped execution.

## 1. Goal

Resolve the CH-dominated fast spectrum without changing the model, initial
equilibrium, P2 space, mesh, quadrature, mobility or saved delta=1e-3 bump.
The sequence and numerical thresholds were fixed in [STEP3A3_DESIGN.md](STEP3A3_DESIGN.md)
before PDE runs. No first interval or t=0 power is omitted. Historical reports
and step3/step3a1/step3a2 artifacts are unchanged.

## 2. STEP3A.2 unresolved temporal problem

Initial excess energy was 5.542888467657825e-8 J/m and D_CH(0) was
0.11400313905309394 W/m. Thus tau_E=4.862048987156724e-7 s; old dt/tau_E
were 20.56746, 10.28373, 5.1418651.
The old nonzero transient FAIL is not rewritten. Its roundoff-scale initial
semidiscrete and last-interval BE identities motivate a time-resolution
investigation, not a free-energy or weak-form change.

## 3. Linearized CH derivation

M q_dot + mobility K r = 0, M r = H q, hence

`L = -mobility M^(-1) K M^(-1) H`.

Only sparse matrices and solves are used at full size. H is obtained with UFL
as the second derivative of the existing discrete scalar free energy.
This analysis does not replace the production CHNS solver.

## 4. Discrete matrices M/K/H

There are 37153 phase DOFs, degree-12 quadrature, P2 and the
unchanged 96x48 crossed mesh. H includes bulk
`lambda/epsilon*(3*phi_e^2-1) q chi`, the gradient term and the exact wall term.
Current `f_wall=-sigma*cos(theta)*(3phi-phi^3)/4`, so
`f_wall''=3*sigma*cos(theta)*phi/2`; no wall Hessian approximation is used.

| matrix | Frobenius symmetry defect |
| --- | --- |
| M0 | 4.139877e-17 |
| K | 9.4847824e-17 |
| H | 6.8348425e-17 |

## 5. Mass-zero projection

`Pq=q-1_h*(m^T q)/(m^T 1_h)`, with assembled `m_i=integral N_i`.
Projection is permitted only inside analysis; there is no production mass fix.
Twenty deterministic vectors give max relative mass defect
3.058468593e-17 and relative `PLPq-Lq` defect
2.454061418e-16. Both pass 1e-11.

## 6. Finite-difference operator validation

An independent nonlinear weak chemical assembly produces the semidiscrete rate,
without time advance. Central differences remove the O(delta) nonlinear term.
The one-sided limit is reported too, not confused with central accuracy.

| delta | Hessian M-dual relative defect | central rate relative L2 | one-sided rate relative L2 |
| --- | --- | --- | --- |
| 0.001 | 2.0121816e-07 | 3.1105761e-09 | 9.827063e-06 |
| 0.0005 | 5.0304551e-08 | 7.8055234e-10 | 4.9135314e-06 |
| 0.00025 | 1.2576167e-08 | 2.3457755e-10 | 2.4567661e-06 |
| 0.000125 | 3.1443229e-09 | 2.6156681e-10 | 1.2283827e-06 |

Three additional deterministic random FE directions also PASS. The final
central-rate plateau is roundoff; the mandatory 1e-6 operator gate passes.

## 7. CH-dominated early-time justification

Read-only historical early samples give max `(D_visc+D_slip)/D_CH` =
1.503294371e-08
and max `E_kin/initial_excess` =
4.827716787e-10.
The spectral planner analyzes the CH-dominated stiff subsystem, while final
qualification remains a FULL nonlinear CHNS experiment, not performed here.

## 8. Relevant Ritz spectrum

The own M-inner-product, twice-reorthogonalized Arnoldi was seeded with the
historical psi coefficients. Dimensions 32/48/64 were increased through 1024.
All 33 converged significant Ritz values are real and negative;
their range is [-10735011121.188166, -10581814086.703129] s^-1.
Not every approximate Ritz value is a converged eigenmode. The full spectrum
is not claimed resolved from this subset.

Initial saved residuals were terminal Arnoldi recurrence estimates. An
independent actual full-operator audit of all 33 converged significant vectors
passed, max residual 2.696463207e-08.
It replaces sub-roundoff estimates in the red plot markers and converged table
rows, without rewriting the saved primary spectrum. Unconverged table rows
still show labelled recurrence estimates. The current API uses full residuals.

| Ritz mode | Re(lambda) [1/s] | Im | tau [s] | relative residual | old z(2.5e-6) | new z(dt0) | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | -3849.3605 | 0 | 0.00025978341 | 27.895775 | 0.0096234013 | 8.9645006e-09 | False |
| 945 | -9788271.2 | 0 | 1.0216309e-07 | 0.69124243 | 24.470678 | 2.2795205e-05 | False |
| 55 | -1.0581814e+10 | 0 | 9.4501755e-11 | 2.6964632e-08 | 26454.535 | 0.02464323 | True |
| 3 | -1.0735011e+10 | 0 | 9.3153141e-11 | 1.0060026e-14 | 26837.528 | 0.025 | True |

Complete records: `validation_results/step3a3/spectrum/ritz.csv`.
Figure: [spectrum.pdf](validation_results/step3a3/spectrum.pdf),
[spectrum.png](validation_results/step3a3/spectrum.png).

## 9. Krylov self-convergence

The single polynomial reference FAILED: m=768/1024 max relative L2 difference
0.14478141.
This is not a detected physical instability: orthogonality and the Arnoldi
relation remained near machine precision. `rho*t_end` is about 1.07e6.

The installed SLEPc 3.24 MFN alternative was investigated without installing
packages. A tiny nonidentity-M test rejected direct BV/M normalization; the
correct Euclidean-internal adapter passed its dense test. Its full-horizon
ncv=32 run reached t=7.711254e-5 at 2736.52 s and was interrupted at the
one-hour resource guard before completing the horizon. It is NOT qualified.
See the primary [SLEPc MFN documentation](https://slepc.upv.es/release/documentation/manual/mfn.html)
and [SciPy matrix exponential documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html)
for the underlying tooling; numerical evidence here comes from saved runs.

A separately versioned rational M-Krylov reference passed levels 32/48/64.
Fixed positive resolvent times are geometrically spaced from 1/rho to 1e-4 s.
The sparse BE resolvent only generates a basis; reference evolution is
`V exp(t V^T M L V) e1 * ||q0||_M`, NOT a BE time history. Balanced auxiliary
scaling and pinned MUMPS change only linear algebra. No model or gate changed.

| 48/64 comparison | maximum relative difference | threshold |
| --- | --- | --- |
| max_relative_L2_difference | 3.6160886e-07 | 1e-06 |
| max_relative_energy_difference | 4.7148441e-10 | 1e-06 |
| max_relative_power_difference | 1.322432e-09 | 1e-06 |

M-orthogonality defect is 2.049124275e-15.
The rational reference's individual Ritz vectors are not asserted converged;
its exponential action is self-converged in full FEM L2/E2/D2. Fast relevant
rates for the planner come from the independently converged primary Ritz set.

## 10. Semidiscrete linear energy identity

`E2=q^T H q/2`, `D2=mobility*(M^-1 Hq)^T K(M^-1 Hq)`.
The full operator's max relative `dE2/dt+D2` defect is
4.152471876e-15; the qualified reference's
projected identity defect is 4.534196074e-11.
Both pass the 1e-10 identity gate. Integrated reference power is evaluated
independently from its quadratic modal integrand, never replaced with -DeltaE.

## 11. Initial spectral timescales

At the unchanged delta=1e-3, linear E2(0)=5.542887728317465e-08
J/m, D2(0)=0.1140031365946258 W/m and
tau_E_linear=4.862048443480074e-07 s. This agrees with the historical
nonlinear timescale. The fastest converged significantly seeded rate is
rho=10735011121.18817 s^-1, tau=9.315314057069357e-11 s.
Its weight is small, but exceeds the predeclared significance criterion; this
is not an unseeded all-grid spectral-radius estimate. Energy timescale and the
fastest weakly excited timescale are different quantities.

## 12. Low-mode temporal benchmark

An independently refined full-operator mode has lambda=-3259.575709365278
s^-1, tau=0.0003067883949211063 s, residual
1.284023952e-09. Its M norm is 1;
max initial |delta phi|=1e-4. This different input is explicitly the isolated
low-mode benchmark, NOT a replacement for the fixed stiff bump.

| steps | dt | measured amplitude | linear BE amplitude | error to exp(-1) | integral D_CH | DeltaE | budget defect | defect/\|DeltaE\| |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 1.533942e-05 | 0.3768922 | 0.37688948 | 0.0090127592 | 5.0039725e-10 | -4.8790139e-10 | 1.2495862e-11 | 0.02561145 |
| 40 | 7.6697099e-06 | 0.3724334 | 0.37243062 | 0.0045539561 | 4.960003e-10 | -4.8980126e-10 | 6.1990478e-12 | 0.012656251 |
| 80 | 3.8348549e-06 | 0.37016959 | 0.37016679 | 0.0022901466 | 4.9384448e-10 | -4.907565e-10 | 3.0879792e-12 | 0.0062922838 |

Observed orders: [0.9848484713382726, 0.9916804555448524]. Amplitude-halving relative BE
errors: [7.210411702152086e-06, 3.6052641558814713e-06, 1.8027630996897415e-06]. Max mass/domain error over all
five low-mode runs is 5.952028993e-14;
minimum certified resolution is 8.328131076.
All isolated low-mode gates PASS, including integrated D and continuum budget
convergence. This does not qualify full CHNS on the stiff bump.

## 13. Spectral timestep planner

z_target=0.05 gives initial candidate dt0=4.657657028534679e-12 s.
The candidate is proposed entirely from the saved linear input; blocks grow
by at most two, with a 1e-4 initial-D2 content cutoff used only for proposal.
The first proposal failed strict global linear error gates. Uniform halving
passed: dt0=2.328828514267339e-12 s, first-block z=0.025,
single-step BE modal amplification error=0.0002998440692.
Both candidates, including the rejection, are retained. No nonlinear observed
field can modify a frozen schedule, and an unauthorized plan is rejected by
the execution-schedule API.

## 14. Planned schedule

This is an immutable **reduced-linear-verified candidate, NOT an authorized
production plan**. Total 1796 intervals to 1e-4 s; half/quarter plans would use
3592/7184 intervals, identical boundaries. `active z` uses content still above
the proposal cutoff; the complete JSON also records rho_initial*dt for every
block (including already-decayed modes).

| block | t_start | t_end | dt | steps | predicted max L2 error | active stiffness [1/s] | active z |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0 | 5.6823416e-10 | 2.3288285e-12 | 244 | 1.0459007e-06 | 1.0602942e+10 | 0.024692434 |
| 1 | 5.6823416e-10 | 1.0805764e-09 | 4.657657e-12 | 110 | 8.9558581e-07 | 4.7061832e+09 | 0.021919787 |
| 2 | 1.0805764e-09 | 1.9375853e-09 | 9.3153141e-12 | 92 | 1.0719805e-06 | 2.4524652e+09 | 0.022845484 |
| 3 | 1.9375853e-09 | 3.3162518e-09 | 1.8630628e-11 | 74 | 1.6624872e-06 | 1.2345145e+09 | 0.02299978 |
| 4 | 3.3162518e-09 | 7.3404675e-09 | 3.7261256e-11 | 108 | 2.8795559e-06 | 5.8642826e+08 | 0.021851054 |
| 5 | 7.3404675e-09 | 1.4494629e-08 | 7.4522512e-11 | 96 | 4.409232e-06 | 2.9027624e+08 | 0.021632115 |
| 6 | 1.4494629e-08 | 2.9697221e-08 | 1.4904502e-10 | 102 | 7.0678065e-06 | 1.3862001e+08 | 0.020660623 |
| 7 | 2.9697221e-08 | 4.9371165e-08 | 2.9809005e-10 | 66 | 1.0596886e-05 | 66372593 | 0.01978501 |
| 8 | 4.9371165e-08 | 1.0541209e-07 | 5.961801e-10 | 94 | 1.7797786e-05 | 40723941 | 0.024278803 |
| 9 | 1.0541209e-07 | 2.2464811e-07 | 1.1923602e-09 | 100 | 2.6838395e-05 | 19592669 | 0.023361519 |
| 10 | 2.2464811e-07 | 4.106563e-07 | 2.3847204e-09 | 78 | 3.8312059e-05 | 9364625.8 | 0.022332014 |
| 11 | 4.106563e-07 | 8.4944486e-07 | 4.7694408e-09 | 92 | 5.6041129e-05 | 4294645.8 | 0.020483059 |
| 12 | 8.4944486e-07 | 1.1928446e-06 | 9.5388816e-09 | 36 | 7.5789222e-05 | 2103018.3 | 0.020060443 |
| 13 | 1.1928446e-06 | 2.7572212e-06 | 1.9077763e-08 | 82 | 0.00013423904 | 1202400.2 | 0.022939106 |
| 14 | 2.7572212e-06 | 4.5123754e-06 | 3.8155526e-08 | 46 | 0.00016825316 | 579040.62 | 0.0220936 |
| 15 | 4.5123754e-06 | 7.1069512e-06 | 7.6311053e-08 | 34 | 0.00024736116 | 257286.45 | 0.0196338 |
| 16 | 7.1069512e-06 | 1.8095743e-05 | 1.5262211e-07 | 72 | 0.00054612729 | 145542.46 | 0.022212997 |
| 17 | 1.8095743e-05 | 9.9901191e-05 | 3.0524421e-07 | 268 | 0.00092044647 | 59534.556 | 0.018172578 |
| 18 | 9.9901191e-05 | 0.0001 | 4.9404333e-08 | 2 | 0.0007788654 | 13442.393 | 0.00066411244 |

File: [timestep_plan.json](validation_results/step3a3/stiff_bump/planner/timestep_plan.json).
File SHA256: `1ad91330d9d7501699f3bbd2d8551db1d26607867ac3bb649642f9fd61d620a3`.
Semantic plan SHA256: `5300ed2a4ca1112d5812a569f6a83ec764c33074d65a1bab7ed868770acfe817`.

## 15. Linear BE schedule verification

| predicted reduced error | value | threshold |
| --- | --- | --- |
| endpoint_relative_L2_error | 0.00077833888 | 0.001 |
| max_checkpoint_relative_L2_error | 0.00092044647 | 0.001 |
| max_quadratic_energy_relative_error | 0.00082851566 | 0.001 |
| integrated_D2_relative_error | 0.0017172779 | 0.01 |

Reference power integral at delta=1e-3 is 4.326537415911412e-08 J/m;
scheduled reduced BE predicts 4.333967283164525e-08 J/m.
The plan's stored unscaled q=psi powers/integrals must be multiplied by delta².
Sparse full-size BE verification is implemented and tested against a dense
tiny tangent system, but its scientific schedule run was **NOT RUN after cost STOP**.
Reduced verification alone does not authorize production nonlinear execution.

## 16. Full nonlinear stiff-bump level 0

NOT RUN. Historical conservative seconds/accepted interval =
9.318659643, giving
4.648975755 hours for level0 and
32.54283029 hours for all three
levels. Predeclared limits are 3 and 10 hours. They were not increased after
seeing the candidate. The measured blocker is computational cost, not a newly
observed weak-form defect. A separate zero-time audit reloaded exactly the
frozen field and re-ran only the existing consistent weak mu initialization:
D_CH(0)=0.1140031390530939 W/m,
excess energy=5.542888467657825e-08 J/m. The actual
mixed-state phi hash matches the historical initial hash exactly. This audit
has ZERO accepted time intervals and does not bypass the cost STOP.

## 17. Half-plan refinement

NOT RUN: level0 not authorized. The deterministic rule and 3592-interval
schedule exist; no measured full nonlinear refinement claim is made.

## 18. Optional quarter-plan refinement

NOT RUN: level0 not authorized. The deterministic rule and 7184-interval
schedule exist, but cost blocks the series and endpoint/integral order evidence.

## 19. Integrated dissipation convergence

| level | min dt (planned) | planned steps | accepted | integral D_CH | DeltaE | final defect | defect/\|DeltaE\| | end phi difference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 2.3288285e-12 | 1796 | 0 | NOT RUN | NOT RUN | NOT RUN | NOT RUN | NOT RUN |
| 1 | 1.1644143e-12 | 3592 | 0 | NOT RUN | NOT RUN | NOT RUN | NOT RUN | NOT RUN |
| 2 | 5.8220713e-13 | 7184 | 0 | NOT RUN | NOT RUN | NOT RUN | NOT RUN | NOT RUN |

The full-CHNS successive D_CH differences and final <=5% gate are unmeasured.
The low-mode PASS and linear integral are not substituted for these gates.

## 20. Energy closure

Full nonlinear final budget, observed budget order and endpoint convergence
remain unmeasured. For audit only, the historical last-pair Richardson hint is
4.13927765e-08 J/m; it is NOT qualification evidence.
The new linear reference is consistent with the hypothesis of temporal error,
but cannot confirm the nonlinear continuum limit by itself.
Figure: [stiff_startup_energy.pdf](validation_results/step3a3/stiff_startup_energy.pdf)
shows the qualified linear reference and labelled historical data only;
t=0 is visible on a symlog axis, not removed.

## 21. Early BE work decomposition

Physical continuum dissipation is **D_CH + D_visc + D_slip**.
**B_BE is a signed numerical/discrete BE work remainder, diagnostic only**.
It can contain negative quartic/cubic pieces and is never clipped or added to
the physical power curve. The identity is
`continuum_local = weak_work_defect - B_BE + dt/2*(D_old-D_new)`.

Measured isolated low-mode first interval (20-step run):

| quantity | value |
| --- | --- |
| Delta_E | -5.2870195e-11 |
| D_old | 3.707313e-06 |
| D_new | 3.3626237e-06 |
| trapezoid_integral | 5.4224363e-11 |
| BE_endpoint_integral | 5.1580696e-11 |
| weak_work_defect | 5.6535353e-19 |
| B_BE | 1.2894995e-12 |
| trapezoid_gap | 2.6436669e-12 |
| continuum_local_defect | 1.3541681e-12 |
| decomposition_roundoff | 6.2202972e-20 |

Thus even this resolved low mode has a visible signed BE/trapezoid split while
weak work is near roundoff. Every low-mode interval was audited. A pure stiff
z=100 regression also has D_old/D_new=10201 with a large trapezoid error and
roundoff BE work. New full stiff-bump first-20 audits are **NOT RUN**.
Figure: [early_be_work_decomposition.pdf](validation_results/step3a3/early_be_work_decomposition.pdf)
is labelled isolated CH, never presented as full stiff-bump evidence.

## 22. Linear-vs-nonlinear comparison

Available: t=0 nonlinear/linear timescale agreement and measured low-mode BE
decay versus its linear prediction. Unavailable: common-early-time full-CHNS
bump L2/H1/excess-energy/power comparison. No new nonlinear snapshot exists.
Low mode answers the isolated modal time-integrator question; stiff bump would
answer practical resolution of the broad fast spectrum. They remain separate.

## 23. Remaining limitations

The cost STOP requires a separately agreed compute resource or implementation
cost reduction before full sparse verification and full nonlinear series. No
claim is made that the candidate is runtime-optimal. Full CHNS production
schedule execution and its first-20 work audit remain unperformed, and no
small-dt production conditioning conclusion can be inferred without that run.

No spreading, theta sweep, epsilon/mobility/slip study, film, tank, unmatched
density, gravity or water-air validation was attempted. There is no continuum
or full model qualification here.

Tests: 250 pure passed, 14 skipped
(FEniCSx/PETSc/SLEPc unavailable on host; 2 expensive tests deselected by the
existing default configuration); 149
Docker tiny checks passed. Existing negative-resolution smoke emits its
expected warning. `.github/workflows/step3.yml` includes new tiny tests,
never the production 96x48 series. GitHub-hosted Actions were NOT run: no push
or remote workflow dispatch. XML evidence is in `validation_results/step3a3/tests/`.

Provenance: base git `29d40b6cbddb882d5b6248dc2f52df5f8f93b2eb`, pinned image
`sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`. Each completed stage stores
source hash/archive, equilibrium/mesh/input fingerprints; plan stores policy
and Krylov hashes. Early interrupted reference prototypes retain source
archives and explicit terminal status but not complete final provenance.
The qualified rational basis/reference, sparse matrices and exact recovered
input arrays are retained. The failed 1024-vector basis is local-only, with
all numerical summaries and hashes retained in git.

Prepared equilibrium fingerprint: `df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc`.
Psi coefficient SHA256: `831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0`.
Initial phi SHA256: `6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681`.
Both initial hashes match historical STEP3A.2 exactly. The historical bump
had only hashes, not a standalone saved psi array; recovery replayed the
unchanged initializer once, verified BOTH hashes and froze the identical bytes.

## 24. Scientific verdict

**MODEL NOT YET VALIDATED.**

PASS: sparse operator/finite differences/conservation/symmetry, relevant fast
Ritz evidence, rational exponential reference, linear energy identity,
isolated low-mode temporal convergence and reduced-linear candidate accuracy.

STOP: predeclared production cost limit. NOT RUN: full sparse candidate
verification and full nonlinear stiff-bump plan/half/quarter convergence.
Therefore **BASIC NONZERO CH TRANSIENT TEMPORALLY QUALIFIED** is not claimed.
Historical STEP3A.2 FAIL remains unchanged. No physics or thresholds were
relaxed, and no numerical BE remainder was relabelled as physical heat.
