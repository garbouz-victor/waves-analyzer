# STEP3A10 — numerical-identity guard and schedule-aware cost audit

**Administrative numerical guard PASS. All four cost backtests PASS.
Continuation NOT AUTHORIZED: full-L0 forecast 16354.593801770276 s exceeds
the unchanged 14400 s cap. No interval 128 was executed.**

The earlier STEP3A10 source-guard STOP and all STEP3A9 reports/results remain
unchanged historical records. This report covers the separately authorized
administrative resumption. Overall: **MODEL NOT YET VALIDATED**.

## 1. Goal

Separate immutable numerical identity from execution lineage, then test the
mandated exact-event cost model before considering completion of full CHNS L0.
Cost planning is not itself a numerical-physics qualification.

## 2. Why STEP3A9 stopped

The historical conservative whole-step p95 model forecast 15942.5699042599 s
for complete L0 against 14400 s. It remains a valid historical rejection under
that policy. Its number, STOP and pilot level0_authorized=false were not changed.

## 3. Frozen scientific trajectory and new guard

The user authorized a new external guard comparing numerical identity while
recording the current session HEAD separately. Exactly one field is excluded
from numerical comparison: `execution_base_HEAD`. All remaining implementation
keys, stack, policies, design, schedule, source archive and audit proof are
checked. All 72 archived numerical source/runner/design files match byte-for-byte.

| provenance | value |
| --- | --- |
| historical execution base HEAD, unchanged | 5f2d9498e1d0715f57b125ee49b4ddfd533467ae |
| actual current session HEAD at audit | 2e7ff2f80c6ce23eb4a7a72e02ffde79999c70d3 |
| numerical implementation SHA | 64250f222cb8fe4eda98e2abd7ec0f91be3d0045bd65a45e81f098833795b717 |
| numerical identity SHA excluding lineage | 3c5750c879403b190ec02f9b88592d80348034c649cc48007410ca2061ba090e |
| old production source archive SHA | f81b460f8a73f7cddbcd0f77e6427601b207419abb39bfa48d55568138ad9f6e |
| old production policy SHA | 557896045044c8af06efe09744944a0549160051e64cfde5e02b808cbae9935d |

The original guard, wrapper, environment, solver sources and checkpoint 127
are not modified or monkeypatched. The new controller passes the new guard
through the frozen DivergenceTrajectory callback API and would use the same
CHNSPhaseRateBE/runner classes. No numerical reimplementation is introduced.
Its numerical execution branch was NOT RUN because cost authorization failed.
The audit's dirty administrative source is bound by individual file SHA hashes;
a later final commit must not replace the recorded actual execution HEAD.

## 4. Whole-step p95 hypothesis and measured limitation

The proposed improvement was to count expensive scheduled archival events
only where they actually occur. The schedule is indeed sparse. However, the
measured long archival tail is NOT confined to those scheduled events:
22 of 109 ordinary F0_C0_A0 observations have archival_s above 1 s, spanning
6.6963–8.8728 s. This 1 s description is diagnostic, not a new cost class or gate.

Examples: step 7 archival 7.362417 s; step 15 archival 6.961270 s; step 69
archival 8.872828 s. All have field=false, checkpoint=false, audit=false.
The frozen ordinary archive path includes the scalar journal append, JSON,
flush/fsync. The existing aggregate timing cannot identify which operation
or storage/OS effect caused the stalls. No new attribution is claimed.

The ordinary archival p95 is 8.370449625197216 s. Applying the prescribed rule
therefore still assigns 10.46306203149652 s of archival allowance to each ordinary
future step. No outliers were discarded or reclassified from their timing.
Separately conservative core/archive estimators make the new total larger
than the historical whole-step estimate. The old model is not called wrong.

## 5. Exact L0 event schedule

The unchanged L0 descriptor has 1068 intervals; schedule SHA:
`21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31`.
Flags come only from field_indices, checkpoint_indices and audit_indices.
The exact count table was saved and printed BEFORE any cost estimate.

| archive event class | total 1–1068 | observed 1–127 | future 128–1068 |
| --- | --- | --- | --- |
| F0_C0_A0 | 1023 | 109 | 914 |
| F0_C1_A0 | 5 | 0 | 5 |
| F1_C0_A0 | 23 | 6 | 17 |
| F1_C0_A1 | 12 | 8 | 4 |
| F1_C1_A0 | 3 | 3 | 0 |
| F1_C1_A1 | 2 | 1 | 1 |
| total | 1068 | 127 | 941 |

Core classes: normal 118 observed/936 future; audit 9 observed/5 future.

## 6. Core/archive cost decomposition

All observed core_s=total_s-archival_s are valid and nonnegative. No material
negative-accounting error, NaN or inf occurred. The timing journal's hash chain
and exact step mapping were verified. Core classification depends only on the
audit flag; field/checkpoint flags do not split the core cost distribution.

## 7. Frozen cost estimator and component table

Policy SHA `476dae1d71bacd7a25df47ec175cf67000525ee9cc2335860440ab3e0b707d1d`.
The policy and administrative source hashes were frozen before calibration and
before examining any forecast/backtest result. n>=5 uses 1.25*p95 with linear
interpolation; n=2..4 uses 1.50*max; n=1 uses 2.00*sample. An unseen required core
class would fail; none occurred. No factors were tuned after measurement.

| cost component | n observed / rule | estimator [s] | future count | forecast [s] |
| --- | --- | --- | --- | --- |
| normal core | 118 / 1.25*p95 | 6.148459278986138 | 936 | 5754.957885131025 |
| audit core | 9 / 1.25*p95 | 6.144880204748915 | 5 | 30.724401023744576 |
| ordinary archive F0_C0_A0 | 109 / 1.25*p95 | 10.46306203149652 | 914 | 9563.23869678782 |
| checkpoint archive F0_C1_A0 | 0 / calibration | 2.1083054009977786 | 5 | 10.541527004988893 |
| field archive F1_C0_A0 | 6 / 1.25*p95 | 10.87800357562628 | 17 | 184.92606078564677 |
| field+audit F1_C0_A1 | 8 / 1.25*p95 | 10.074433957501697 | 4 | 40.297735830006786 |
| field+checkpoint F1_C1_A0 | 3 / 1.50*max | 2.2079616989940405 | 0 | 0 scheduled contribution |
| field+checkpoint+audit F1_C1_A1 | 1 / 2.00*sample | 3.8585550020070514 | 1 | 3.8585550020070514 |
| initialization | 1.25*max completed primary init | 10.365156259995274 | 1 | 10.365156259995274 |
| session overhead | 1.50*max overhead | 1.993266609035345 | 1 | 1.993266609035345 |

Future interval subtotal 15588.544861565239 s; with initialization/overhead:
15600.90328443427 s. No sparse event is treated as free.

## 8. Sparse-event calibration

Two unseen-in-required-training cases each received exactly 5 replays.
F0_C1_A0 is unseen in the complete prefix and needed for five future checkpoints.
F1_C1_A1 is observed at 122 for the main forecast, but unseen in the relevant
backtest training subsets; its calibration is used ONLY for those subsets.
The final main forecast retains the required 2*observed-at-122 estimate.

| repetition | F0_C1_A0 [s] | F1_C1_A1 [s] |
| --- | --- | --- |
| 1 | 1.3993573259940604 | 1.7469120529931388 |
| 2 | 1.3061254610001924 | 1.798956790997181 |
| 3 | 1.405536933998519 | 1.740951210995263 |
| 4 | 1.3754725210019387 | 1.659145254001487 |
| 5 | 1.3490661870018812 | 1.7881641720014159 |
| 1.50*max allowance | 2.1083054009977786 | 2.6984351864957716 |

The actual frozen archive, atomic_json and Journal paths were used with copied
checkpoint 127 arrays, same /work bind filesystem (device 123), compression and
fsync unchanged. Field payload hashes reproduce the existing field 127 layout.
Checkpoint replay includes the numerical guard, journal-prefix calculation,
checkpoint NPZ/metadata and LATEST write.

Implementation detail: the full frozen runner uses expanded_audit=False for
the isolated-only residual file path; its actual full audit is embedded in
scalar rows. Consequently no fictitious extra residual JSON file was added or
existing full audit removed. The audit flag remains a distinct cost class.

Payload descriptions, sizes, hashes and all timings were recorded before
deleting only the temporary calibration copies. The physical trajectory was
not copied into a new trajectory or modified. Entire calibration-function
wall 17.77648990399757 s, including setup/final record I/O/cleanup, is conservatively
included in total planning. The final pre-cleanup calibration.json elapsed
snapshot is 17.417705533000117 s; these clock boundaries are explicit.

## 9. Rolling backtests

| window | predicted interval wall [s] | actual interval wall [s] | predicted/actual | pass |
| --- | --- | --- | --- | --- |
| 51–64, train 1–50 | 207.4671706251288 | 117.15722545600875 | 1.77084400742343 | PASS |
| 65–96, train 1–64 | 525.4084758797981 | 194.72125168195635 | 2.6982595445614863 | PASS |
| 97–127, train 1–96 | 513.6498803194721 | 150.77114419701684 | 3.406818214818822 | PASS |
| complete prefix 1–127, train 1–50 | 1874.6223399532687 | 734.8168520499676 | 2.551142280860203 | PASS |

Every validation window is upper-bounded without tolerance. Passing backtests
does NOT override a failed full-L0 budget gate.

## 10. Schedule-aware forecast

Observed session 0: wall 281.41992658700474, init 8.137124692002544,
intervals 272.16723071498564, overhead 1.115571180016559 s.
Session 1: wall 472.2705907490017, init 8.29212500799622,
intervals 462.64962133498193, overhead 1.3288444060235634 s.
Both overheads are finite and positive. Actual prefix wall 753.6905173360065 s
is used, not a fitted replacement.

| quantity | old forecast [s] | schedule-aware forecast [s] | actual complete [s] |
| --- | --- | --- | --- |
| complete L0 | 15942.5699042599 | 16354.593801770276 | N/A, incomplete |
| remaining 941 intervals/session | 15188.879386923894 | 15600.90328443427 | N/A, not run |
| total scientific allowance usage | 16181.252679854915 | 16611.053067269288 | N/A, no full run |

Historical scientific spent 992.3732929310208 s includes audit, preparations,
moderate/tiny/restart and primary prefix exactly once. New total prediction adds
calibration 17.77648990399757 s and future continuation 15600.90328443427 s.
The component table explains the difference; no new timing distribution was
substituted for an observed class.

## 11. Authorization decision

| mandatory condition | measured | unchanged cap | pass |
| --- | --- | --- | --- |
| complete L0 forecast | 16354.593801770276 s | 14400 s | FAIL |
| total scientific forecast | 16611.053067269288 s | 19800 s | PASS |
| four backtests | all upper-bound actual | no underprediction | PASS |

**authorized=false**, without slack or manual override.
Decision SHA `d2c8d93fc1e338a197556b837e4778e54b85f1d55b4ab2155c370e53a72a1b2a`.
Authorization-file SHA `f4c41892840e92f74c7117b58a680365a0fbaf3c43c8feebd0b8a997d36576ff`.
The false authorization still binds policy, decision, manifest, schedule,
old numerical policy/archive and all administrative files.

## 12. Prefix provenance before continuation

The existing manifest SHA remains
`ec8421ff9f7f28c84229294fcf17cffddbae464e398efbefe76aac64c5cf4127`.
Scalar history SHA `5916da0860220fb0fa785c208897cbb6cab062fc465cada5ee4a1e0cac7f6a91`;
timing history SHA `e7fb4fcce29ca6c6db6e5c1d96c789d06d06e5187dd592c0c281191ded8a95fe`.
Both hash chains and all historical file bindings PASS before and after audit.
Latest accepted checkpoint 127 semantic SHA
`bc229cc210b61b6b5ef374f9dd33ed69b88f540065652bc3b415c306e7222d21` is unchanged.
No row 0–127 was rewritten or reserialized.

## 13. Continuation execution

NOT RUN. The actual controller was invoked against authorized=false and refused
before importing frozen PDE classes. No session directory, solver or interval 128
was created by that refusal check. Starting accepted state is still 127 at
t=6.148107277665776e-10 s; there was no actual continuation start.
Controller SHA `f6709a7d87063567eb0a84ef702344d16ed0c92ada6babce9b587f26c6e1b215`.
The separate numerical guard PASS did not bypass cost authorization.

## 14. Actual cost

Primary full-path wall remains historical 753.6905173360065 s. No additional
primary PDE time was spent. Historical total plus conservatively charged
administrative calibration is 1010.1497828350184 s. General tests/report rendering
are separate administrative work, not disguised primary trajectory wall.
Actual complete-L0 wall and forecast/actual ratio are N/A.

## 15. L0 completeness

**INCOMPLETE 127/1068; no COMPLETE.json; no FAILED.json.** No new rejected
primary interval. Zero rejections over the existing 127 are inherited evidence,
not proof of a zero-rejection complete 1068-step trajectory. No partial completion
or disconnected probe is substituted for the missing full horizon.

## 16. Full physical gates

No new physical state was produced. Complete-L0 physical extrema are NOT RUN.
The unchanged prefix has max Newton 1, max mass/domain 1.9737298215558338e-14,
min certified cells 8.328131075834218 and every recorded physical gate passing.

## 17. Divergence over full L0

Complete-L0 divergence statistics NOT RUN. Inherited prefix, including t=0:
strong min/p50/p95/max 0 / 3.263664111387802e-8 / 5.615297465930906e-8 /
6.044641155438896e-8; max weak continuity 2.2840662037059463e-23;
max projection fraction 4.979369080447153e-16. These values were not hidden,
promoted to full-horizon results, or compared to the retired old strong cap.

## 18. Full energy/work

| requested complete-L0 quantity | result |
| --- | --- |
| Delta E_total | NOT RUN |
| integrated D_CH / D_visc / D_slip / D_total | NOT RUN |
| final budget defect / budget divided by actual energy change | NOT RUN |
| sum B_BE / max weak-work / decomposition / work identity | NOT RUN |

No energy correction, missing t=0 power sample, or B_BE-as-heat interpretation.

## 19. Full-vs-isolated phi coupling

NOT RUN at full horizon. No extrapolation of tiny prefix coupling to 1e-4 s.

## 20. D_CH coupling

Full-horizon integrated-D coupling NOT RUN. Existing trapezoidal power
accounting is unchanged; it is not replaced by minus energy change.

## 21. Free-energy coupling

Final full F_CH and F coupling NOT RUN; precision floor unchanged.

## 22. Kinetic/hydrodynamic coupling

Full-horizon maxima and cumulative hydro dissipation NOT RUN. No full-horizon
physical curves were produced. New PDFs contain cost observations only.

## 23. Transfer vs isolated temporal uncertainty

| observable | full coupling | inherited isolated L0–L2 uncertainty | threshold | pass |
| --- | --- | --- | --- | --- |
| max common phi | NOT RUN | 3.691188008762449e-8 | coupling<=gap and initial-relative<=1e-3 | NOT RUN |
| final integrated D_CH | NOT RUN | 7.033312454099817e-11 J/m | coupling<=gap | NOT RUN |
| final F_CH | NOT RUN | 6.0652038946784614e-12 J/m | coupling<=max(gap,1e-13) | NOT RUN |
| Ekin/excess | NOT RUN | N/A | <=1e-4 | NOT RUN |
| significant hydroD/DCH | NOT RUN | N/A | <=1e-4 | NOT RUN |
| cumulative hydroD/excess | NOT RUN | N/A | <=1e-4 | NOT RUN |

Isolated gap values are inherited committed references, not a new comparison.
The required post-COMPLETE reference re-verification/transfer analysis did not
run because full L0 was not authorized. Other frozen D/F significance thresholds
also remain unchanged and unevaluated at the missing horizon.

## 24. Combined error proxies

Full-horizon combined phi and integrated-D proxies NOT RUN. They would be
practical conservative proxies, not rigorous continuum error bounds.

## 25. Tests, figures and remaining limitations

New combined STEP3A10 preflight: 30 PASS on host and 30 PASS pinned. Before freeze,
older host NumPy rejected the newer `method=` argument; the API call was made
compatible without changing linear-percentile mathematics. Initial failure XML
is retained; no measured cost data existed at that correction.

Final complete host suite: **323 passed, 46 skipped, 2 deselected**, 50.81 s.
Pinned STEP3–STEP3A10 suite: **254 passed**, 82.45 s, one expected warning from
an intentionally underresolved negative fixture. New tests cover lineage versus
numerical identity, exact classes/counts, every estimator branch, unseen classes,
sparse-event inclusion, initialization/overhead, negative accounting, both caps,
backtest failure, hash binding and refusal before PDE imports. The workflow
adds only tiny/pure tests, never the production continuation.

No remote CI for the unpushed source was observed. No push. Cost PDFs:
cost_event_classes.pdf, cost_forecast_backtest.pdf, forecast_vs_actual.pdf.
The last figure omits actual complete-L0 because it is unavailable.

The remaining planning issue includes ordinary journal archival latency, not
only scheduled field/checkpoint events. No I/O/compression/fsync/schedule change
or alternative latency model was implemented after seeing the result. Such work
would need a separate explicitly authorized iteration. No moving contact-line,
sensitivities, coupled temporal levels, film or tank work was started.

## 26. Verdict

Administrative identity separation: **QUALIFIED BY READ-ONLY HASH CHECKS**.
Cost backtest conservatism: **PASS ON ALL FOUR PRESCRIBED WINDOWS**.

**SCHEDULE-AWARE COST AUDIT DID NOT AUTHORIZE COMPLETION WITHIN THE FROZEN
14400s L0 CAP.** The total 19800 s cap passes, but cannot override the L0 cap.

**FULL-CHNS COUPLING TRANSFER REMAINS INCOMPLETE DUE TO COST GATE.**

**MODEL NOT YET VALIDATED.**
