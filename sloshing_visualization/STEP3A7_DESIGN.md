# STEP 3A.7 — fixed multi-step isolated-CH protocol

Execution base: dfe272049123ecbd9d200af8ff070ae4eb4b08d3. Overall:
**MODEL NOT YET VALIDATED**. This design, machine policy, complete runner,
checkpoint/restart and analysis sources are frozen before the scientific
pilot. No source changes are permitted within a scientific series.

## Inherited evidence and immutable inputs

STEP3A6 qualified moderate algebraic/root accuracy and one exact tiny BE
step, NOT a trajectory. Read-only historical inputs: STEP3A2 prepared
equilibrium, STEP3A3 saved perturbation/operator, STEP3A4 optimized and
full-sparse-verified plan, STEP3A6 proof. No spectral/optimizer rerun.

- initial phi SHA: 6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681
- psi SHA: 831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0
- equilibrium fingerprint: df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc
- parent plan SHA: 5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8

All actual initialized coefficients must match. D_CH(0)=0.11400313905309394
within the inherited relative 1e-10 check; certified cells reproduce
8.328131075834218 within relative 1e-12 and remain >=8. Physical fingerprint,
mesh, P2/P2, quadrature 12, delta=1e-3, matched density, zero body force,
free/wall energy, mobility/epsilon/sigma/theta/slip are unchanged. No field
correction, clipping, smoothing or artificial energy/dissipation accounting.

## Schedules and guesses

Expand every parent block with its original `t_start+k*dt` endpoints; the
committed end equals that expression. L0 retains each original solver dt
bitwise and all 1068 intervals. L1/L2 split EACH L0 physical interval [a,b]
at a+k*(b-a)/factor, explicitly assigning the last child endpoint to b.
Child solver dt is parent solver dt/factor (factors 2,4); no new block grid
or optimizer. This is an exact mathematical subdivision with the unavoidable
floating-point clock/nominal-dt discrepancy recorded separately. Existing
Diagnostics integrates powers using differences of accepted physical times.
No clock rounding is used to adapt a solver dt. Parent endpoints, no gaps,
strict monotonicity and final t=1e-4 are checked. Expected counts 1068/2136/4272.
Every expanded schedule has its own semantic SHA and parent/derivation version.

Every interval is phase-rate BE. Controls throughout: atol=1e-10, rtol=1e-9,
stol=0, max_it=30, newtonls/bt, preonly/LU/MUMPS; actual getters checked each
attempt. First guess: existing unprojected semidiscrete rate. Subsequently:
previous accepted rate, NEVER rescaled at a dt change. No retry, alternative
guess, tolerances, max_it, line search, absolute BE or BDF2 path.

## Physical checks and persistent accounting

Same existing material checker and mandatory EVERY accepted-state checks:
finite arrays/metrics; mass relative to initial AND target / domain <=1e-10;
exact/certified P2 cells >=8; exactly two wall crossings; energy growth
<=1e-9 J/m. Every interval uses existing BEWorkDiagnostics, with absolute
weak-work and decomposition limits 1e-12 J/m. These values are taken from
committed STEP3A4 policy JSON, not relaxed. The initial power is retained.

Persistent accepted Diagnostics.initial/previous/cumulative is transactionally
copied for candidates. Cached diagnostic forms refer to fixed private physical
Functions. Only after SNES/material/physical/work success are physical and
diagnostic histories published. Failed attempts retain accepted fields,
history, cumulative work and last rate unchanged; rejected fields are separate.

Primary physical I_D is the independent trapezoid integral including t=0.
Existing energy normalizations and symmetric diagnostic are retained. Separate
diagnostics accumulate endpoint dt*D_new, signed B_BE, weak defect and stable
work-minus-remainder energy changes. B_BE is NOT physical heat. Literal
absolute residual is NEVER an acceptance veto.

## Sparse observations, restart and durable storage

Scalars, work, solver/rate/timing statistics and worst certified-cell witness:
one append-only, hash-chained JSONL per level, flushed and fsynced each step.
Full phi/mu/actual retained-rate fields: common L0 endpoints at powers of two,
all block ends, final, plus 25/26/27/28/50 and the first post-boundary five
states/pilot endpoint. The same times are exact child endpoints on L1/L2.
At t=0 the archived rate is labelled semidiscrete initial rate, NOT an
accepted BE unknown. No temporal field interpolation.

Restart checkpoints: initial, step 25/50, first block end, pilot end,
every 200 accepted intervals, and final; orderly prefix stops also checkpoint.
They contain state/old/older, current phi/mu, actual last rate, maps, time/index,
physical/algorithm/source/schedule/stack fingerprints, diagnostic histories,
work sums and hash/byte count of the committed scalar prefix. Atomic directory
publication and atomic LATEST/COMPLETE markers distinguish durable checkpoints
from partial writes. COMPLETE requires the exact final interval, coherent
history and all common field samples. It never means temporal qualification.

On crash recovery, any history/field/checkpoint suffix newer than the selected
durable checkpoint is preserved in a labelled orphan/recovery directory before
replay with identical controls. A scientific FAIL marker prohibits resume.
Fingerprint mismatch, concurrent writer or corrupt hashes prohibit reuse.
Completed levels cannot be continued. Level0 pilot is stored directly in the
level0 dataset and continued, not recomputed; pilot/ contains summaries only.

Expanded read-only residual audit: first interval, powers-of-two common
indices, first block end/first following interval, and final. Retain UFL/
sparse rate, absolute evaluators and chemical/Riesz diagnostics as in STEP3A6.
Forms, mesh-only maps and wall topology are reused; mathematical diagnostics
and frequency of hard checks are not weakened.

## Pilot and restart authorization

First run exactly 50 L0 intervals. Then a diagnostic fork from checkpoint 25
executes 26..28, compared with primary saved states. It does not replace any
primary state. Frozen restart limits: relative M-L2 phi/mu/rate <=1e-10;
cumulative D relative <=1e-10 (denominator max(abs(reference),1e-30)); absolute
F/budget difference <=1e-12. Tiny separate-process restart tests additionally
cross a dt change. Timing differences are excluded from physical-history
equivalence, but retained in each run's own record.

Extend pilot to max(50,first block last index+5), derived from plan: 127 here.
No new initial state. First post-boundary convergence/physical gates required.
Iteration explosion is predeclared as post-boundary Newton count greater
than max(5,2*max(previous five counts)); it stops authorization for investigation.
The optional alternate-guess boundary fork is NOT elected. No further forks.

## Cost limits and forecast (wall time, not CPU)

Source of limits: validation_results/step3a4/policy/policy.json:
L0 <=7200 s, whole isolated series <=28800 s. These limits stay fixed.
The old code's forecast mixed absolute low-mode timing with sparse work
sampling; it cannot describe this phase-rate/every-step-work runner. The
machine policy fixes limits but no percentile rule. For this iteration,
freeze the requested measured-only conservative estimator:

    projected seconds/remaining interval = 1.25*max(mean, p95)

using all available accepted phase-rate interval total times. Include actual
spent wall time, initialization and diagnostic fork costs; never double-count
the reused prefix. Record SNES, diagnostics, resolution, work and I/O separately.
The scalar journal contains pre-commit category timings; a separate append-only
timing journal measures the complete interval including scalar/field/checkpoint
I/O. Only its total is used by the cost estimator. Its own final timing-record
write is included in session wall time, not recursively in that same record.
Reassess after the boundary pilot, complete L0, and complete L1. L0 authorization
also requires the entire projected three-level series to fit; before each
subsequent level require spent+remaining forecast <=28800. Actual runtime guards
remain active during execution. Cost stops are incomplete evidence, not PASS.

## Three-level hard convergence criteria

All L0/L1/L2 COMPLETE markers, matching source/physics/initial hashes, accepted
counts, every-step gates, pilot/boundary/restart checks are prerequisites.
Each refined level starts independently at the original t=0 state.

Retain the strict STEP3A2 energy_convergence definitions (no new roundoff
exception or relaxed trend): absolute successive I_D gaps decrease; last
gap/abs(I2)<=.05; absolute final budget defects strictly decrease; endpoint
M-L2 phi gaps strictly decrease; final F gaps decrease; nonzero energy loss
and positive physical integral. Report also gap/abs(finer integral) for both
pairs. Zero/nonfinite denominators do not pass. Finest final abs(B)/abs(DeltaF)
<=.05 WITHOUT using the shared significance floor to hide this nonzero gate,
as in STEP3A2. Existing shared energy validation is also retained. No hard
5% requirement on each local first interval. Observed orders and max common-time
field/power/integral differences are diagnostic, not fitted requirements.

No third complete level => no formal temporal qualification; no Richardson
extrapolation substitute. No boundary-guess adaptation or lower solver accuracy.

## Frozen implementation and limited execution

Implement and run pure/tiny FEM tests BEFORE source freeze. Archive full
multiphase sources, CLI/analysis runner source hash, design and policy SHA,
execution base HEAD, actual image/PETSc/DOLFINx/MUMPS/precision/rank. Verify
unchanged fingerprints at initialization, checkpoints and continuation. A code
bug after scientific pilot stops this series; results with different hashes
cannot be combined. Historical data and execution provenance are read-only.

Sequence: tests/freeze -> 50-prefix -> fork restart -> boundary prefix -> cost
authorization -> resume L0 -> cost -> fresh L1 -> cost -> fresh L2 -> analysis.
Reports/PDFs/summary are generated from measured histories, including honest
NOT RUN/incomplete statuses. No source editing during the scientific sequence.

Maximum success: **NONLINEAR ISOLATED CH STIFF TRANSIENT TEMPORALLY QUALIFIED
ON THE CERTIFIED MESH; FULL CHNS TRANSFER AND REMAINING STEP 3A PHYSICAL
BENCHMARKS STILL REQUIRED.** Overall **MODEL NOT YET VALIDATED**.
No full CHNS, coupling probe, moving contact, film, tank, gravity/water-air,
unequal density, new spectrum/optimizer, or STEP3B work in this iteration.
