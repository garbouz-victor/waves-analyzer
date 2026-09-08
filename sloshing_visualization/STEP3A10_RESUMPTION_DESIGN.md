# STEP3A10 — authorized administrative resumption

This adds to, and does not rewrite, the first-hard-gate STOP recorded in
STEP3A10_DESIGN/REPORT and the STEP3A9 historical cost STOP. The user explicitly
authorized a separate numerical-identity guard and current-session HEAD record.

## Administrative identity, before any cost analysis

The new external guard reproduces every old source_guard check, except that
`execution_base_HEAD` is separated from numerical implementation comparison.
No other field is ignored. The old guard and every file under multiphase remain
byte-identical. The wrapper environment is read, never set or spoofed.
Historical identity in checkpoint127 is immutable trajectory-origin provenance;
the actual new session HEAD is recorded separately, not merged over old fields.

The old execute() necessarily calls the old HEAD guard internally. The new
controller therefore uses the SAME frozen CHNSPhaseRateBE and
DivergenceTrajectory classes through their existing public callback interface,
with only external session administration and the newly authorized guard.
It does not subclass/patch the numerical runner, solver, or diagnostics. It
duplicates the old execute's initialization and wall accounting administratively;
all PDE steps, controls, transactions, archives, and COMPLETE creation remain
the existing frozen classes' responsibility. This is the limited consequence
of the user's new authorization, not a numerical source revision.

The new numerical guard and prefix127 integrity must PASS before cost analysis.
The original prefix manifest is reused unchanged. Existing STEP3A10 STOP and
summary records are historical; new resumption results are stored separately.

## Frozen mandatory cost rules

Archive classes are (field,checkpoint,audit), from the actual schedule sets;
core classes depend ONLY on audit. Count all steps1..1068, observed1..127,
future128..1068. Freeze the policy before estimating times/calibration/backtests.
Core is total_s-archival_s; negative values below -1e-9s fail, smaller roundoff
is explicitly recorded and set to zero. Nonfinite times fail.

For either observed class: n>=5 uses1.25*p95, with numpy linear-percentile
interpolation; n2..4 uses1.50*max; n1 uses2.00*sample. Missing required core
class fails. Missing required archive class gets exactly5 read-only archival
replays, then1.50*max. No factor tuning or selecting repetitions by outcome.

Calibration uses copied checkpoint127 arrays and frozen archive/atomic_json/
Journal code on the SAME /work bind filesystem, in a temporary directory
under step3a10/cost_audit. It retains compression/fsync and replays scalar,
field, checkpoint+LATEST, and residual-JSON operations for the required class.
Checkpoint replay includes guard and journal-prefix accounting. Existing full
trajectory has no separate isolated residual file output (expanded_audit=False);
the full audit results are in scalar rows, so replay follows that actual code.
Unseen classes required by any backtest use this identical mechanism, not
future observed samples. Payload descriptions/hashes/sizes/times are saved
before deleting temporary calibration files. No PDE solve is allowed.

Initialization allowance=1.25*max completed primary initialization_JIT_s.
Session overhead=wall-init-sum(interval total_s), from exact session endpoints;
overhead<-1s fails, otherwise1.50*max(overhead,0). There is one forecast session.
Future cost is the sum of exact class counts times their class estimators, plus
initialization and session overhead. Actual primary prefix wall is added once.
Total scientific spent includes audit and all preliminary sessions without
double-counting the primary prefix. Read-only cost calibration administrative
wall is reported separately and conservatively added to total planning.

Backtests are fixed: train1..50/predict51..64;1..64/65..96;1..96/97..127;
and train1..50/predict1..127. Actual interval totals must not exceed predictions
in ANY test. Missing core, failed accounting, underprediction or either cap
failure means STOP; no post-hoc estimator modification or continuation.

Full-L0 cap14400s and total scientific cap19800s are unchanged, without slack.
Decision/policy/manifest/schedule/numerical policy/archive and administrative
source SHA bindings form the separate authorization before interval128.

## Conditional execution and analysis

Before continuation: recheck authorization, numerical guard, exact127 prefix,
LATEST and absence of FAILED/COMPLETE. Current HEAD and controller SHA are
recorded in the new session. Reuse the immutable inherited trajectory identity
without calling it the current execution HEAD. Resume exactly128..1068, once.
Actual remaining wall=min(14400-used full,19800-used total), minus initialization;
the frozen run_until wall enforcement is active. No second session after cap
or scientific failure. No retries, guess/dt/physics/diagnostic/logging changes.

Only after the frozen runner creates a valid COMPLETE: use existing read-only
comparison functions against reverified isolated L0/L2, same scalar/common field
times, all unchanged significance/dominance/energy gates. Record full extrema,
energy/work, combined practical proxies, actual versus forecast and figures.
An incomplete prefix never qualifies transfer. MODEL NOT YET VALIDATED always.

New pure tests cover numerical identity vs lineage, policy estimators, exact
event counting, unseen-class handling, cap/backtest failures and authorization
bindings. No production trajectory is placed in CI.
