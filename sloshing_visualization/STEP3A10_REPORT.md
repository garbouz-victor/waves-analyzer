# STEP3A10 — first-hard-gate checkpoint

**STOP BEFORE COST ANALYSIS.** The actual unchanged STEP3A9 `source_guard()`
fails on execution-base HEAD identity. No event-cost forecast, calibration,
backtest, authorization, or new PDE interval was performed.

**MODEL NOT YET VALIDATED.** This is neither evidence against the proposed
schedule-aware cost estimator nor a newly failed physical/PDE step.

## 1. Goal

Qualify the mandatory schedule-aware cost model within the unchanged 14400s
full-L0 and 19800s total caps, then conditionally complete the same trajectory.
The user requires an exact source_guard PASS before any cost analysis.

## 2. Why STEP3A9 stopped

The historical conservative whole-step p95 model projected 15942.5699042599s
for full L0, exceeding 14400s. The historical STOP, pilot authorization=false,
reports, cost policy and forecast remain unchanged. No claim that the old
model was wrong is made.

## 3. Frozen scientific trajectory and current STOP

Execution repository HEAD was `c718347d81925faa9a6a8184f4dd4c2d0298aaae`.
The unchanged Docker wrapper obtains that value with host `git rev-parse HEAD`
and exports it as STEP3_GIT_COMMIT. `implementation()` includes this variable
as execution_base_HEAD, and `source_guard()` compares the entire implementation
dictionary with the frozen record.

| guard operand | frozen | current at preflight |
| --- | --- | --- |
| execution_base_HEAD | 5f2d9498e1d0715f57b125ee49b4ddfd533467ae | c718347d81925faa9a6a8184f4dd4c2d0298aaae |
| multiphase source SHA | 64250f222cb8fe4eda98e2abd7ec0f91be3d0045bd65a45e81f098833795b717 | identical |
| runner SHA | 00163e66320f118de853b3de07cee5ab453759dc8976ff9c6aefd5e03d96e839 | identical |
| production archive SHA | f81b460f8a73f7cddbcd0f77e6427601b207419abb39bfa48d55568138ad9f6e | identical |
| production policy SHA | 557896045044c8af06efe09744944a0549160051e64cfde5e02b808cbae9935d | identical |
| design, inherited policy, audit proof, schedule, stack | frozen values | identical |

All archived numerical source/runner/design members are independently checked
byte-for-byte, and all match. All six solver-core hashes match. The **only**
implementation difference is execution_base_HEAD. This is a provenance
identity incompatibility across the already committed checkpoint, not evidence
of changed numerical equations or stack. The original guard nonetheless
returns `ValueError: STOP: source/policy/stack/schedule changed after freeze`.
Adding new STEP3A10 files outside the frozen source tree did not change any
numerical source hash or this outcome.

No environment variable was overridden, no guard monkeypatched, no frozen
record edited, and no Git reset performed to manufacture a PASS. Resolving
this requires explicit direction on inherited numerical identity versus the
new session's actual HEAD; the current restrictions do not authorize that.

## 4. Why whole-step p95 is conservative

The requested hypothesis concerns multiplying sparse-event whole-step costs
by every remaining interval. It is a proposed planning improvement, not a
measured STEP3A10 result. It was not tested after the first hard gate failed.

## 5. Exact L0 event schedule

L0 count and descriptor integrity PASS: 1068 intervals; SHA
`21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31`.
The frozen descriptor is bound in inherited_prefix_manifest.json. Event-class
enumeration and observed/future cost-class counts are NOT RUN, because they
belong to the downstream cost analysis. No event schedule was modified.

## 6. Core/archive decomposition

NOT RUN. Timing journal integrity was checked, without computing new timing
statistics or decomposing total_s minus archival_s.

## 7. Cost estimator

The prescribed 1.25*p95 / 1.50*max / 2.00*sample rules and class definitions
are recorded as an inactive downstream plan in STEP3A10_DESIGN.md. No factor
or budget was tuned, and no active cost policy/decision was created.

| cost component | observed estimator | future count | forecast |
| --- | --- | --- | --- |
| normal core | NOT RUN | NOT RUN | NOT RUN |
| audit core | NOT RUN | NOT RUN | NOT RUN |
| ordinary archive | NOT RUN | NOT RUN | NOT RUN |
| field archive | NOT RUN | NOT RUN | NOT RUN |
| checkpoint archive | NOT RUN | NOT RUN | NOT RUN |
| audit/combined archive classes | NOT RUN | NOT RUN | NOT RUN |
| initialization | NOT RUN | one planned session, not authorized | NOT RUN |
| session overhead | NOT RUN | one planned session, not authorized | NOT RUN |

## 8. Sparse-event calibration

NOT RUN. No calibration copies, temporary I/O payloads, compression changes,
or trajectory writes occurred.

## 9. Rolling backtests

| validation window | predicted | actual window cost | ratio | pass |
| --- | --- | --- | --- | --- |
| 51–64, train1–50 | NOT RUN | NOT ANALYZED | N/A | NOT RUN |
| 65–96, train1–64 | NOT RUN | NOT ANALYZED | N/A | NOT RUN |
| 97–127, train1–96 | NOT RUN | NOT ANALYZED | N/A | NOT RUN |
| 1–127, train1–50 | NOT RUN | NOT ANALYZED | N/A | NOT RUN |

No backtest failed or passed; no factor retuning was attempted.

## 10. Schedule-aware forecast

| quantity | historical whole-step forecast [s] | schedule-aware [s] | actual complete [s] |
| --- | --- | --- | --- |
| complete L0 | 15942.5699042599 | NOT RUN | N/A |
| remaining L0 | 15188.879386923894 | NOT RUN | N/A |
| total scientific wall | 16181.252679854915 | NOT RUN | N/A |

Historical numbers are quoted from STEP3A9, not recalculated or changed.
Caps remain 14400s and 19800s. Whether the new model fits is unknown.

## 11. Authorization decision

Continuation **NOT AUTHORIZED**: first source_guard prerequisite failed.
No continuation_authorization.json, cost decision, or executable continuation
controller was created. The preflight cannot call execute(), full_L0(), or
authorize a continuation even if a future preflight itself passes.

## 12. Prefix provenance before continuation

Read-only diagnostic checks following guard failure found:

| prefix property | result |
| --- | --- |
| accepted endpoint | 127/1068 |
| physical time [s] | 6.148107277665776e-10 |
| scalar history records | 128, steps0–127; valid chain |
| timing records | 127, steps1–127; valid chain |
| scalar history SHA | 5916da0860220fb0fa785c208897cbb6cab062fc465cada5ee4a1e0cac7f6a91 |
| timing history SHA | e7fb4fcce29ca6c6db6e5c1d96c789d06d06e5187dd592c0c281191ded8a95fe |
| checkpoint127 semantic SHA | bc229cc210b61b6b5ef374f9dd33ed69b88f540065652bc3b415c306e7222d21 |
| checkpoint127 NPZ SHA | d60fd633db9a79d86cff8b23e863eb8c21db184aef0198ef39d10bea1341f189 |
| LATEST | correct accepted checkpoint127, valid SHA |
| physical configuration / inputs / mesh | match inherited frozen identity |
| fields and checkpoints | all saved array/metadata hashes verified |
| history clock / controls / recorded hard gates | PASS |
| primary sessions | both closed complete |
| FAILED.json / COMPLETE.json | neither exists |

Manifest SHA:
`ec8421ff9f7f28c84229294fcf17cffddbae464e398efbefe76aac64c5cf4127`.
It includes file sizes and hashes, session records, policies, schedule descriptor,
archive identity and retained rate. It is immutable read-only provenance, not
cost or execution authorization. No row0–127 was reserialized or copied into
a new physical trajectory.

## 13. Continuation execution

NOT RUN. Interval128 was not attempted. Final accepted endpoint remains127.
No PDE retry, new solver, dt/control/guess change, or changed state publication.
No new primary session was started.

## 14. Actual cost

Actual complete-L0 wall and forecast/actual ratio: N/A, L0 incomplete.
No new scientific trajectory wall was spent. Administrative read-only
preflight and unit tests are not reported as scientific PDE execution.

## 15. L0 completeness

INCOMPLETE, 127/1068. No COMPLETE marker. No manually fabricated marker.
Zero new rejected attempts; the inherited prefix has zero rejected intervals.

## 16. Full physical gates

Full-horizon extrema NOT RUN. Inherited prefix evidence remains unchanged:
max Newton1, max mass/domain1.9737298215558338e-14, min certified cells
8.328131075834218. These are historical prefix results, not full-L0 claims.

## 17. Divergence over full L0

NOT RUN over complete L0. Historical prefix max weak continuity is
2.2840662037059463e-23 and max f_Q4.979369080447153e-16. Strong divergence
min/p50/p95/max, including initial u=0, was
0 / 3.263664111387802e-8 / 5.615297465930906e-8 / 6.044641155438896e-8.
No old empirical hard cap was reintroduced.

## 18. Full energy/work

| full-horizon quantity | result |
| --- | --- |
| Delta E_total | NOT RUN |
| integrated D_CH / D_visc / D_slip / D_total | NOT RUN |
| final budget / budget divided by actual energy change | NOT RUN |
| sum B_BE / max weak-work / decomposition / work identity | NOT RUN |

B_BE remains a numerical BE remainder, not physical heat. Prefix energy
closure must not be substituted for final full-horizon closure.

## 19. Full-vs-isolated phi coupling

NOT RUN at full horizon; no transfer gate evaluated from prefix data.

## 20. D_CH coupling

NOT RUN at full horizon. Integrated physical power was not replaced by energy
loss; no history samples were removed.

## 21. Free-energy coupling

NOT RUN at full horizon; the energy precision floor remains unchanged.

## 22. Kinetic/hydrodynamic coupling

Full-horizon Ekin/excess, significant hydroD/DCH and cumulative hydroD/excess
are NOT RUN. No continuation curves or new physical figures were generated.

## 23. Transfer vs isolated temporal uncertainty

| observable | full coupling | inherited isolated L0–L2 gap | frozen threshold | pass |
| --- | --- | --- | --- | --- |
| max common phi L2 | NOT RUN | 3.691188008762449e-8 | coupling <= gap; initial-scale <=1e-3 | NOT RUN |
| final integrated D_CH | NOT RUN | 7.033312454099817e-11 J/m | coupling <= gap | NOT RUN |
| final F_CH | NOT RUN | 6.0652038946784614e-12 J/m | coupling <= max(gap,1e-13) | NOT RUN |
| Ekin/excess | NOT RUN | N/A | <=1e-4 | NOT RUN |
| hydroD/DCH, significant | NOT RUN | N/A | <=1e-4 | NOT RUN |
| cumulative hydroD/excess | NOT RUN | N/A | <=1e-4 | NOT RUN |

Reference numbers are explicitly inherited from the committed STEP3A9 audit.
The requested post-COMPLETE isolated-reference re-verification and transfer
analysis were not run because L0 did not complete.

## 24. Combined error proxies

Full-horizon combined phi and integrated-D proxies: NOT RUN. No extrapolation
of negligible prefix coupling to t=1e-4 is claimed.

## 25. Remaining limitations and tests

To proceed, the administrative execution-HEAD identity must be resolved with
explicit user direction while keeping numerical source, physics and historical
provenance intact. No fix is implemented under this iteration's restrictions.

New first-gate/provenance tests: local host **7 passed**; same pure tests in
the pinned container **7 passed**. These test fail-closed status and byte/size
bindings, not the unimplemented schedule-aware estimator or a PDE continuation.
The actual pinned source_guard itself **FAILS**, as required to report honestly.
No full FEM regression or new production experiment was run. Remote CI for
these unpushed changes was not observed. No workflow change was needed.

An initial invocation of the new administrative CLI tried querying Git inside
/work, which does not mount .git; it failed without writing artifacts or calling
a solver. The CLI was corrected to read the unchanged wrapper's existing
STEP3_GIT_COMMIT environment value (never set/override it), then recorded the
actual guard failure and verified prefix. This was not a scientific retry or
a changed cost model. The earlier direct source_guard call already showed the
same execution-HEAD mismatch before any new file was created.

All moving-contact, sensitivity, falling-film and tank work remains out of scope.
No push has been performed.

## 26. Scientific verdict

Cost: **SCHEDULE-AWARE COST AUDIT NOT RUN — BLOCKED BY FROZEN SOURCE-GUARD
EXECUTION-HEAD MISMATCH; CONTINUATION NOT AUTHORIZED.**

Transfer: **FULL-CHNS TRANSFER REMAINS INCOMPLETE; STEP3A10 STOPPED AT
SOURCE-GUARD PREFLIGHT.** The inherited scientific PASS results stand.

**MODEL NOT YET VALIDATED.**
