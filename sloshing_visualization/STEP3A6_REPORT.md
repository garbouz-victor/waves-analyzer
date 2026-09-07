# STEP 3A.6 — matched Newton accuracy and target tiny-step proof

## 1. Outcome and scope

Both questions of this iteration have measured answers:

1. **The moderate mismatch disappears after matching nonlinear root accuracy.**
   The unchanged phi/mu/D comparison limits of 1e-10 pass at both P1 and P2.
2. **The exact target rate-form step converges at the original production
   tolerances**, with either semidiscrete or zero initial rate. A fresh
   absolute negative control reproduces the historical line-search failure.

**TARGET TINY-dt ISOLATED CH BE STEP SOLVED;
ABSOLUTE-STATE CANCELLATION HYPOTHESIS SUPPORTED.**

Overall verdict: **MODEL NOT YET VALIDATED**.
This is not a pilot, temporal-convergence series, full-CHNS transfer or
moving-contact validation. No such run was authorized or performed.

The [design](STEP3A6_DESIGN.md) and
[machine policy](validation_results/step3a6/policy/policy.json) were frozen
before new scientific PDE calculations. Historical STEP3A4/5 results,
their reports and their FAIL verdicts were not changed.

## 2. Base, immutable inputs and provenance

Execution HEAD was `7ed5609cb9f13012f750b3f7fce7c98454386b4c`, with the new
implementation in the worktree. Every scientific stage archives its actual
source and records that HEAD separately from historical input provenance.
The measured implementation source SHA256 is
`b03c1b0b50462467793d1237362572f7d572f3b30d300e9d9f02ff4ba01767f2`.
Prior artifacts are not attributed to this implementation.

Pinned stack: DOLFINx 0.10.0, PETSc 3.24.0, MUMPS 5.8.1, UFL 2025.2.0,
real float64; image
`sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`.
Serial mesh: 18,432 crossed triangles, 37,153 DOFs per scalar P2 field;
quadrature degree 12. The original IsolatedCH, CHNSSolver, BDF2, continuum
equations, material/free/wall energies, physical configuration and mesh were
not modified. No clipping, smoothing or mass/energy correction was used.

The actual initialized arrays and frozen plan were checked, not merely their
metadata labels:

| input | SHA256 / fingerprint |
|---|---|
| initial phi | `6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681` |
| psi | `831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0` |
| prepared equilibrium | `df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc` |
| optimized 1068-step plan | `5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8` |
| mesh | `93dbd2cf76148bcd50afef9e353c55f92c15789ad0e303dbb97435a4672b9b15` |
| initial consistent mu | `e711921b7399300a33c9875d3d4e1f59025a072dc0c9715515c58f8ee03bcbad` |

Physical fingerprint `d42900144d7d75c67e4455351ac98a0c41083623ed06de5c2987c1d35a485650`
is separate from validation-only nonlinear controls. The embedded observer
config is not the one-step experiment schedule: the actual interval/dt is
recorded explicitly. Its untouched default SNES fields are likewise not the
active moderate overrides; the actual getter-based controls are a separate
object. Source, input, initial physical guess, unknown and
reconstructed-state hashes are retained in each run directory.

## 3. Algebraic equivalence, Newton accuracy and temporal accuracy

The original phase residual is

    R_abs(s) = ((phi_new-phi_old)/dt,s) + mobility*(grad(mu),grad(s)).

The alternate residual and reconstruction are

    R_rate(s) = (r_phi,s) + mobility*(grad(mu),grad(s)),
    phi_new = phi_old + dt*r_phi.

Both use exactly `(mu,chi)-F'_h(phi_new)[chi]=0`. At corresponding trial
states in exact arithmetic, **the residual rows are identical**. The change
of unknown transforms columns: `J_rate=J_abs*diag(dt,I_mu)`. It does not
automatically rescale residual rows. Small FE tests verify assembled
residual/Jacobian equivalence (relative 1e-12) and a finite-difference action.

Primary initial guesses are deliberately independent: absolute starts at
phi_old; semidiscrete-rate starts at phi_old+dt*r_semidiscrete. This explains
different initial residuals without invoking residual-row scaling. A raw
mixed Euclidean residual is not the relative error of either physical field.

All moderate solves independently cover **[0,1e-6 s]** from the same phi_old
and mu_old. No previous endpoint became another solve's old state; no solution
was copied between formulations. The moderate experiment checks algebra and
Newton accuracy, not BE temporal accuracy. The tiny experiment proves one
step is solvable; it does not establish a resolved trajectory.

## 4. Retained unknowns and transactional acceptance

The phase-rate implementation now retains independent copies of old phi,
the actual mixed phase-rate coefficients, new mu, rounded physical phi,
dt*rate, physical-state increment, recovered-rate diagnostic, mixed/scalar
maps, Function and PETSc vectors. Every array has a SHA256 and an archive
roundtrip check. The actual rate is **never recovered by subtraction for
storage or used interchangeably with the recovered-rate diagnostic**.

Installed DOLFINx `NonlinearProblem.solve` copies Function -> PETSc x before
SNES and x -> Function on normal return. We additionally synchronize using
the installed public `assign` API and verify serial coefficient equality;
we do not assume shared memory. On failure, the last residual Function,
internal SNES iterate and last accepted physical state have separate archive
names. The fresh negative-control Function and SNES iterate happened to be
identical; this is measured, not an API assumption.

Acceptance is transactional: SNES success -> synchronization -> private
reconstruction -> existing material checker -> finite/mass/resolution/
topology/energy/work gates -> publication of physical fields and histories.
Diagnostic cumulative state is private until the candidate passes. Archive
failure also rolls back the one-step validation transaction. Tests include
material rejection after successful SNES with a nonempty accepted history;
state/old/older, time/counters/current dt/phase, last rate, snapshot, history
and cumulative diagnostics remain unchanged.

The original absolute solver source is untouched. Its original problem is
solved through a validation-only candidate/publication wrapper.

## 5. TABLE A — moderate accuracy

Controls P1: atol=rtol=1e-12; P2: atol=rtol=3e-13. All retain stol=0,
max_it=30, newtonls/bt, preonly/LU/MUMPS. Actual getters before and after each
solve confirmed those controls. Every new solve converged with reason +2
in 2 iterations. The final residual blocks below are unscaled algebraic norms.
Differences are consistent mass-norm differences of fields, not differences
of norms; mu's mean was not removed.

| formulation | guess | atol / rtol | iterations | phase residual | chemical residual | phi rel difference | mu rel difference | D_CH rel difference |
|---|---|---|---:|---:|---:|---:|---:|---:|
| P0 absolute, historical | phi_old | 1e-10 / 1e-9 | 1 | 1.14e-13* | 1.77e-11* | reference | reference | reference |
| P0 rate, historical | semidiscrete | 1e-10 / 1e-9 | 2 | not retained | not retained | 1.976e-11 | 1.6725e-8 FAIL | 1.5873e-9 FAIL |
| P1 absolute | phi_old | 1e-12 / 1e-12 | 2 | 1.05979e-13 | 1.10651e-16 | reference | reference | reference |
| P1 rate | semidiscrete | 1e-12 / 1e-12 | 2 | 1.00812e-14 | 3.54954e-15 | 5.11954e-15 | 2.79788e-12 | 1.78746e-12 |
| P2 absolute | phi_old | 3e-13 / 3e-13 | 2 | 1.05979e-13 | 1.10651e-16 | reference | reference | reference |
| P2 rate | semidiscrete | 3e-13 / 3e-13 | 2 | 1.00812e-14 | 3.54954e-15 | 5.11954e-15 | 2.79788e-12 | 1.78746e-12 |
| P2 rate | zero | 3e-13 / 3e-13 | 2 | 9.91768e-15 | 4.03996e-18 | 3.85679e-17 | 2.24661e-15 | 1.08802e-14 |

*P0 absolute blocks are the historical postmortem reassembly. P0 rate's true
unknown was not archived; its recovered-rate residual is not relabelled here
as a retained-rate residual. P0 is historical evidence, not a new solve.

P1/P2 primary |F difference|=4.16334e-17 J/m and |mass difference|=1.11022e-16 m²,
against unchanged 1e-12 limits. P1 and P2 endpoints are coefficientwise
identical within each primary formulation. Tightening absolute P0 -> P1
changes mu by 1.67271e-8 relative and D by 1.58906e-9; rate P0 -> P1 is
unchanged. The old mismatch was therefore limited by the absolute solver's
earlier Newton stopping, not mathematical non-equivalence.

All five new moderate candidates pass distinct physical groups: finite and
material admissible, two wall crossings, certified cells >=8, maximum
mass/domain 3.39235e-15, no forbidden energy growth, and absolute work/
decomposition gates. They are not qualified by agreement alone.

## 6. TABLE B — remaining Newton-error estimates

One sparse diagnostic `J e=-R` at each P2 root; no corrected state was
published. In the rate representation `e_phi=dt*e_rate`. For D, the estimate
uses `2*mobility*mu^T K e_mu` plus its separately recorded quadratic term.

| formulation | level / guess | estimated phi rel error | estimated mu rel error | estimated D relative sensitivity |
|---|---|---:|---:|---:|
| absolute | P2 / phi_old | 3.26886e-17 | 2.00256e-15 | 2.10963e-14 |
| phase rate | P2 / semidiscrete | 5.11959e-15 | 2.79793e-12 | 1.79310e-12 |
| phase rate | P2 / zero | 3.06649e-18 | 8.50352e-16 | 5.91099e-16 |

All satisfy the frozen quarter-threshold engineering targets: relative
2.5e-11 for phi/mu/D, absolute 2.5e-13 for F/mass. Largest estimated absolute
F and mass errors are 3.34625e-20 J/m and 1.86772e-19 m². Correction solves
have positive KSP reason 4; relative linear residuals are respectively
5.44e-15, 5.14e-11, 4.08e-14. The corrections are finite and small. These
are local Newton estimates, **not rigorous a posteriori bounds**.

The absolute chemical algebraic residual decreased from historical 1.77e-11
to 1.11e-16, about 1.6e5-fold. Its current remaining mu-error estimate is
2.00e-15, far below the measured P0 -> P1 change of 1.67e-8. No separate
P0 correction solve was performed and no strict P0 error bound is asserted.

P2 zero-rate and absolute have identical physical initial-guess hashes and
identical initial residual 0.784388568715275. Their final comparison passes.
P1/P2 constancy is expected when the same two Newton iterations are taken;
no artificial convergence order versus solver tolerance is claimed.

## 7. Linear sparse tiny-step comparator

Only after moderate qualification, the historical balanced sparse linear BE
update was compared with

    [ M0             mobility*K ] [rate] = [0   ]
    [-dt*H           M0         ] [mu  ]   [H*q0],
    q_new = q0 + dt*rate.

dt is exactly 4.657657028534679e-12 s. Historical q0 is projected psi
**without delta**; nonlinear delta=1e-3 is applied exactly once. The old
analysis update includes its existing roundoff mass projection; the new
linear rate reconstruction is uncorrected. No nonlinear state is projected.

Endpoint relative L2 difference 2.63009e-16 <=1e-10;
increment relative L2 difference 1.95285e-11 <=1e-6;
mass/domain 1.06183e-18 <=1e-11. Both endpoints, the true linear rate,
chemical field and both increment evaluations are archived with hashes.
KSP reason is 4; sparse system relative residual 9.06e-11.
No spectrum, optimizer or full-plan verification was rerun.

## 8. TABLE C — exact target nonlinear tiny-step

All three attempts use dt=**4.657657028534679e-12 s** and the SAME fresh input,
with actual atol=1e-10, rtol=1e-9, stol=0, max_it=30, newtonls/bt,
preonly/LU/MUMPS. There was no tolerance leakage, timestep change or retry.

| formulation | initial guess | dt [s] | SNES reason | initial / final norm | iterations | accepted | SNES / complete experiment wall [s] |
|---|---|---:|---|---|---:|---|---:|
| rate primary | semidiscrete | 4.657657028534679e-12 | +2 FNORM_ABS | 1.04793e-8 / 1.01695e-14 | 1 | yes | 1.766 / 9.848 |
| rate control | zero | 4.657657028534679e-12 | +2 FNORM_ABS | 7.84389e-1 / 1.01817e-14 | 1 | yes | 1.964 / 11.361 |
| absolute negative control | phi_old | 4.657657028534679e-12 | −6 DIVERGED_LINE_SEARCH | 7.84389e-1 / 1.11091e-8 | 11 | **no** | 20.933 / 26.896 |

Primary effective target is 1e-10 (rtol*r0=1.04793e-17). Zero-rate and
absolute targets are 7.84389e-10 (same initial physical guess/residual).
All requested tolerances are identical despite different initial norms.

The primary uses one full bt step (lambda=1), one direct KSP iteration,
KSP reason 4, MUMPS INFOG(1)=0. Complete iteration and raw line-search
traces are saved. Its Jacobian infinity norm is 10.6667, versus about
6.33665e6 for the absolute attempt; these are **not condition numbers**.
Raw Python lambda fields are null because this pinned binding has no getter;
the raw PETSc monitor supplies the actual accepted/trial-step information.

## 9. TABLE D — one accepted state, different evaluators

These are read-only reassemblies from the primary **retained-rate archive**.
The sparse and UFL assembly paths are independent. Riesz norms use the same
consistent FEM mass matrix: sqrt(R^T M0^-1 R).

| residual evaluator | algebraic norm | L2 Riesz norm |
|---|---:|---:|
| retained rate UFL phase | 1.01694808e-14 | 3.67696817e-12 |
| retained rate sparse phase | 5.11999608e-14 | 1.54695968e-11 |
| original literal absolute UFL phase | 1.28472512e-8 | 3.36699019e-6 |
| coefficient-difference sparse phase | 9.12879915e-9 | 2.28385070e-6 |
| chemical on old+dt*retained rate | 2.81981297e-18 | 9.02081465e-16 |
| chemical on reconstructed physical Function | 5.85960477e-17 | 2.11416043e-14 |

Literal/retained phase norm contrast is **1.26331e6**. The literal value is
also above the unchanged requested SNES target. This old evaluator is the
object of the diagnosis, not a veto on rate acceptance. Conversely, the small
rate residual alone did not establish acceptance: physical, reconstruction,
chemical, archive and independent linear-update checks also passed.

## 10. Reconstruction and ULP audit

Primary retained rate L2 norm =279.03236355834787; max|rate|=15394.503715985304 s^-1;
dt*max|rate|=7.170231843356218e-8. Initial semidiscrete rate L2=283.78413112450124,
max=15786.770344367076 s^-1. These are finite, meaningful physical-rate scales.

The retained increment L2 is 1.2996370493161826e-9. Rounding the physical
coefficients gives increment discrepancy L2=1.0637393399574509e-17,
max=5.5501304319101613e-17, relative=8.184895471525287e-9. Bitwise
reconstruction succeeds, but is not used as a substitute for these error
measurements. An 80-digit Decimal audit of five deterministic worst
coefficients confirms approximately half-ULP reconstruction rounding.
Maximum coefficient ULP/dt=4.767302606540087e-5 s^-1; the mass matrix converts
such coefficient-level uncertainty into the observed algebraic floor.

Retained rate SHA:
`e807cc05ba38a95e98b330f1daaae57284dd5569c1281b2b2d531130450e77d7`.
Rounded new phi SHA:
`9d62fcd07eb1a7d21f3993f29a7ac657331c811160a94d8a62bec545cedf6680`.
Rate, recovered rate and increments are separate arrays in
[state.npz](validation_results/step3a6/target_tiny/rate_semidiscrete/state.npz).

## 11. TABLE E — tiny-step physics and work

Primary accepted state; units are J/m for energy/work and W/m for powers.

| quantity | measured value |
|---|---:|
| mass change / domain | 1.3877787807814457e-15 |
| certified transition cells | 8.328131075834218 |
| D_old | 0.11400313905309394 |
| D_new | 0.11187268676021843 |
| F_old | 0.06564278268454628 |
| F_new | 0.06564278268402256 |
| DeltaF: total-energy subtraction | −5.23719956291302e-13 |
| DeltaF: existing integrand-difference evaluator | −5.235196577654442e-13 |
| DeltaF: existing energy work minus B_BE | −5.235196623222865e-13 |
| physical trapezoid integral | 5.260260638377246e-13 |
| BE endpoint integral | 5.210646057897898e-13 |
| B_BE, numerical remainder only | 2.4550171073465866e-15 |
| weak-work defect | −3.942515010918783e-20 |
| continuum local defect (existing integrand DeltaF) | 2.5064060722803685e-15 |
| trapezoid gap | 4.961458047934762e-15 |
| decomposition roundoff | 4.5568423025713825e-21 |
| weak defect / interval activity, diagnostic only | 7.494904305987054e-8 |
| reconstruction discrepancy L2 | 1.0637393399574509e-17 |
| nonlinear / sparse-linear increment relative difference | 9.991322148389651e-6 |

All existing physical gates pass, including material admissibility, finite
fields, two wall crossings, mass <=1e-10, certified cells >=8, growth <=1e-9
and work/decomposition <=1e-12. D_CH(t=0) is included. Phi DOF extrema are
−0.9790806988061016 and 1.0156953397485375; they were not clipped. The matched
material model remains admissible. Mu mean=0.17992010393728922 and
std=0.0016423668674258069 Pa.

The total-energy subtraction loses relative accuracy on this tiny interval.
All three DeltaF evaluations remain separately labelled. No energy correction
or silent replacement was made. Physical dissipation is D_CH here; velocity,
viscous and slip terms are zero **by isolated-CH construction**, not measured
full-CHNS coupling zeros. B_BE is a signed numerical work remainder and was
never included in physical heat.

For contrast, the correct moderate one-step root still has continuum local
defect about 5.38206e-8 J/m despite weak-work defect about 1e-20. Its coarse
power quadrature is not temporally qualified. Matching Newton accuracy and
qualifying continuum transient accuracy are different tasks.

## 12. Independent linear and rejected-state comparisons

The nonlinear increment compared is **dt*retained rate**, not subtraction
of rounded physical states. The linear comparison uses delta exactly once.
Its relative error 9.99132e-6 is below the unchanged 1e-3 nonlinear comparison
limit; absolute L2 difference=1.29851e-14, physical endpoint relative
difference=1.60959e-14. The retained-linear-increment alternative gives
9.991322125908875e-6. Chemical-update relative discrepancy is 1.41931e-6,
diagnostic only. This agreement is consistent with small nonlinear
corrections; linear reference is not a replacement for the nonlinear solve.

Against the historical STEP3A4 **REJECTED** Function, at matching dt and
initial hashes: phi difference L2=1.59216e-17, maximum coefficient difference
2.22045e-16, mu difference L2=2.49611e-14, D_new difference=−8.98032e-14 W/m,
total-DeltaF difference=1.38778e-17 J/m. The new accepted solution is near
the rejected iterate, but was independently solved, not copied or fitted.

## 13. Matched-guess control and root-cause conclusion

Both rate guesses pass the tiny solver/physical/linear checks in one
iteration. Their physical phi arrays are bitwise equal; mu relative
difference=2.83759e-15 and D relative difference=3.33067e-15.
Zero-rate and absolute have the same initial physical guess and actual
controls. The absolute control reproduces reason −6, 11 iterations and
1.1109110367278936e-8 residual, exactly the historical floor. Its archived
last Function is bitwise equal to the historical rejected mixed Function;
the accepted old physical field remains unchanged.

Together, matched guesses, identical inputs/controls, independent physical
and linear checks, and same-state residual contrast strongly support the
absolute-state cancellation hypothesis. Initial-guess improvement is not
the sole explanation. This remains a double-precision, fixed-case numerical
diagnosis, not a theorem about every mesh/solver configuration.

## 14. Cost, tests and reproducibility

Measured scientific stage wall times total **168.682 s** against the fixed
1800 s budget. These include initialization/JIT, diagnostics, archive I/O
and three diagnostic Newton solves. Small orchestration/JSON aggregation
overhead is not included in that stage sum; there was ample budget margin.
No pilot cost or accepted-trajectory throughput qualification is inferred
from the successful one-step timings.

| scientific stage | total wall [s] |
|---|---:|
| five independent moderate solves + their audits | 111.297 |
| linear sparse pair | 4.932 |
| primary rate tiny | 9.848 |
| zero-rate tiny | 11.361 |
| rejected absolute control | 26.896 |
| target read-only proof arithmetic | 4.347 |

Each run's timing.json separates application categories and inclusive PETSc
events; those overlapping event times must not be summed twice. Full raw
Newton/KSP/MUMPS/line-search traces are retained.

Local tests (not remote CI):

- ordinary suite: **277 passed, 24 skipped, 2 deselected**, 33.05 s;
- pinned STEP3/3A1–3A6 suite: **186 passed**, 81.51 s;
- pre-scientific compact pinned suite: **20 passed**, 2.91 s.

The host skips FEM dependencies; they were exercised in the pinned suite.
The Docker warning is from an existing intentionally under-resolved failure
test, not the production scientific mesh. STEP3 workflow now includes
test_step3a6_*.py; ordinary CI already collects the new pure tests. **Remote
CI was not run/observed for these changes; no push was performed.**

Development failures remain visible in tests/: the first coarse CI fixture
at dt=1e-6 did not converge at strict absolute tolerances; that toy's
moderate interval was set to 1e-3. Scientific dt=1e-6 and the target dt were
never changed. The first full-Docker command did not expand host wildcards
and ran no tests; the final command expanded them inside the container.
Neither event was a retry of a scientific experiment.

CLI: scripts/step3a6_benchmarks.py --stage freeze / moderate-accuracy /
linear-tiny / target-tiny / analyze / all. Scientific output directories
are append-only: occupied experiment paths are rejected. `all` follows only
the bounded prerequisite chain and cannot run pilot/series/full CHNS.
`analyze` reads existing evidence without a PDE solve. Histories, arrays,
controls and source snapshots allow read-only reassembly of this result.

## 15. Explicit NOT RUN and final scientific status

NOT RUN: 50-step pilot; 1068-step level0; half/quarter temporal refinement;
full CHNS rate solver; coupling probe; moving-contact benchmarks; tank,
film drainage, gravity, water-air or unequal-density validation. Historical
P0 moderate was not rerun. No extra tolerance levels, polish, or alternate
dt/line-search attempts were added.

| status | value |
|---|---|
| moderate_accuracy | qualified |
| target_tiny_step | accepted |
| primary nonlinear_solver_status | converged |
| primary physical_checks | passed |
| absolute negative-control execution_status | failed; diagnostic, never accepted |
| cancellation_hypothesis | supported |
| pilot_authorized | false |
| full_series_authorized | false |
| full_CHNS_transfer | not_qualified |

**TARGET TINY-dt ISOLATED CH BE STEP SOLVED;
ABSOLUTE-STATE CANCELLATION HYPOTHESIS SUPPORTED.**

**MODEL NOT YET VALIDATED.** Multi-step nonlinear temporal qualification,
full-CHNS transfer and the remaining STEP3A physical benchmarks are still
required. The iteration stops here.
