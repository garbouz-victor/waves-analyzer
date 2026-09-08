# STEP3A10 — schedule-aware cost qualification

Execution base: `c718347d81925faa9a6a8184f4dd4c2d0298aaae`.
Scope: administrative cost audit and, only after all prerequisites pass,
continuation of the unchanged STEP3A9 trajectory from accepted checkpoint127.
Overall scientific verdict remains MODEL NOT YET VALIDATED.

## Gate 0 precedes all cost analysis

The actual, unmodified STEP3A9 `source_guard()` must pass in the unchanged
pinned Docker wrapper environment. Numerical source equality excluding a
provenance field is not a substitute for this gate. The first direct call
failed before this document was created; no timing estimators, calibration,
event-count table, backtests, or continuation have been executed.

After a failure, only read-only integrity diagnosis and new STEP3A10 STOP
artifacts are permitted. In particular, this iteration must not override
`STEP3_GIT_COMMIT`, replace `source_guard`, edit an old policy/archive, reset
Git HEAD, or reinterpret the old rejected cost decision as authorization.

`scripts/step3a10_preflight.py` records the actual guard outcome, independently
compares its operands, and verifies the prefix/hash chains/accepted archives
without constructing a solver. It writes only new STEP3A10 files and exits
nonzero on guard failure. Manifest creation after failure is diagnostic,
not permission to proceed. Existing sources and historical data are read-only.

## Required downstream plan — NOT EXECUTED OR QUALIFIED

The user's mandated model is retained without modifications: exact schedule
flags define archive classes `(field, checkpoint, audit)` and core classes by
audit only. Core is total minus archival, with negative-roundoff allowance
1e-9s only. Each observed class uses 1.25*p95 for n>=5, 1.50*max for n=2..4,
and 2.00*sample for n=1. Unseen future core fails; unseen archive requires
exactly five same-filesystem, unchanged-format, read-only calibration replays
and 1.50*max. Nothing is calibrated unless it is an unseen required class.

Initialization is 1.25*max measured primary initialization; session overhead
is 1.50*max(nonnegative measured residual session overhead), with material
negative accounting errors rejected. Actual historical wall is used once.
Future event counts multiply their own estimates; sparse I/O is never dropped.

Rolling backtests are 1..50 -> 51..64, 1..64 -> 65..96, and 1..96 -> 97..127.
The fourth test predicts 1..127 using only 1..50 training. Every prediction
must upper-bound measured interval cost. No factor tuning is permitted.

Full L0 cap stays 14400s; total cap stays 19800s. Both forecasts must pass.
Policy, prefix, decision and authorization hashes must be bound before any
interval128 call to the existing `execute(..., resume=True, preliminary=False)`.
Actual caps and all scientific gates remain active. Only the frozen runner
may write COMPLETE; coupling analysis requires a valid complete 1068-step L0.

This downstream plan has no active cost-policy/decision/authorization artifact:
Gate 0 must first be resolved with explicit user direction. No future source
or provenance change is authorized by this document.
