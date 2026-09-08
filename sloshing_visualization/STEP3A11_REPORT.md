# STEP3A11 — Journal / fsync latency attribution

**Historical 7–9 s stalls were not reproduced.** All 480 primary burst appends
and all 60 paced appends were below 1 s. The frozen five-stall minimum was not
met, so neither fsync nor a particular storage path is attributed as root cause.
Interval128 was not run; accepted=127 and the old cost rejection are unchanged.

## 1. Goal

Attribute ordinary journal latency without a PDE solve, cost reauthorization or storage change.

## 2. Inherited frozen scientific state

Prefix 127/1068 at 6.148107277665776e-10 s. Independent guard and all 72 frozen archive members PASS. The historical execution HEAD is preserved; this administrative execution base is f78f90ce2c690ab444de5224c3c395ba2f8a4213.

## 3. Why STEP3A10 cost audit failed

Frozen forecast 16354.593801770276 s exceeded 14400 s; authorization=false remains unchanged. Ordinary archival tails were not confined to sparse events.

## 4. Exact ordinary archival code path

RateTrajectory.advance times scalar Journal.append only on F0_C0_A0. Copy/previous hash → record hash → open(ab) → JSON encode → write → flush → fsync → close → memory publication. Timing journal append is outside both archival_s and total_s. Field/checkpoint NPZ and expanded-audit files are absent inside ordinary archival_s.

## 5. Journal durability semantics

fsync is intentional in the restart model. Atomic_json fsyncs the temporary file then os.replace; archive fsyncs NPZ, publishes metadata, then renames the directory. Neither fsyncs the parent directory. Atomic rename is not by itself a guarantee of persistence after host/VM crash. No durability semantics were changed.

## 6. Attribution policy and predeclared rules

See policy.json and STEP3A11_DESIGN.md. Primary operation decision uses bind burst ONLY: >=5 stalls, median fsync share>=0.90, non-fsync p95<=0.10 s, no other phase median share>=0.20. All repetitions/factors frozen before measurements; no pooling, outlier deletion or retuning. Probe overhead is included conservatively in non-fsync and total append, not the syscall timer.

Frozen policy SHA:
`36037021a673fcf3605b794662bf4424032f4e9385117abebf1cef92ce697444`.

## 7. Environment and mount topology

| context | filesystem | N | stall fraction | fsync p50 [s] | p95 [s] | max [s] |
| --- | --- | --- | --- | --- | --- | --- |
| historical ordinary archival | historical bind | 109 | 0.201834862385 | not phase-timed | not phase-timed | not phase-timed |
| bind | virtiofs | 120 | 0 | 0.00124021049851 | 0.00295129435253 | 0.00478452300013 |
| host | ext4 | 120 | 0 | 0.0010048490949 | 0.00193598791957 | 0.00734421890229 |
| container_tmp | overlay | 120 | 0 | 0.00193777900131 | 0.00829023484621 | 0.00985661300365 |
| docker_volume | ext4 | 120 | 0 | 0.00196577599854 | 0.00820025544745 | 0.0136152019986 |

Docker context: desktop-linux; image sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8. Exact mount sources/options/device IDs, Python versions and complete mountinfo are archived. Filesystem claims use observations, not wrapper comments.

## 8. Historical ordinary stalls

N=109, stalls=22, fraction=0.2018348623853211; statistics: `{"max": 8.872827902996505, "min": 0.003095804997428786, "p50": 0.004499394999584183, "p75": 0.006869280005048495, "p90": 7.714183269398927, "p95": 8.370449625197216, "p99": 8.821156899038469}`. All stall events follow.

| step | physical time [s] | archival [s] | total [s] | core [s] | byte offset | line bytes |
| --- | --- | --- | --- | --- | --- | --- |
| 7 | 3.26035991997e-11 | 7.362416968 | 11.668617522 | 4.30620055401 | 70331 | 10398 |
| 15 | 6.9864855428e-11 | 6.961270124 | 11.17532283 | 4.214052706 | 153626 | 10389 |
| 22 | 1.02468454628e-10 | 7.207489504 | 11.5403255 | 4.33283599601 | 226542 | 10384 |
| 23 | 1.07126111656e-10 | 7.556806013 | 11.016592095 | 3.459786082 | 236926 | 10402 |
| 34 | 1.5836033897e-10 | 7.449055489 | 11.253268477 | 3.80421298801 | 351420 | 10390 |
| 42 | 1.95621595198e-10 | 6.696331211 | 11.041924273 | 4.345593062 | 434550 | 10386 |
| 48 | 2.2356753737e-10 | 7.993548445 | 12.293538289 | 4.299989844 | 496848 | 10386 |
| 54 | 2.51513479541e-10 | 7.961962993 | 13.034927768 | 5.072964775 | 559198 | 10389 |
| 56 | 2.60828793598e-10 | 8.728462356 | 12.896916835 | 4.168454479 | 579994 | 10396 |
| 59 | 2.74801764684e-10 | 8.093529717 | 12.756824112 | 4.663294395 | 611158 | 10393 |
| 61 | 2.84117078741e-10 | 7.701983636 | 12.907507012 | 5.205523376 | 631893 | 10372 |
| 63 | 2.93432392798e-10 | 8.380775998 | 12.968276743 | 4.587500745 | 652645 | 10385 |
| 66 | 3.07405363883e-10 | 7.675354851 | 12.319359564 | 4.64400471299 | 684042 | 10375 |
| 67 | 3.12063020912e-10 | 8.79454390999 | 13.131060142 | 4.336516232 | 694417 | 10396 |
| 69 | 3.21378334969e-10 | 8.872827903 | 13.231810757 | 4.358982854 | 715226 | 10401 |
| 71 | 3.30693649026e-10 | 8.354960066 | 12.484062916 | 4.12910285001 | 735998 | 10383 |
| 80 | 3.72612562283e-10 | 8.823471072 | 13.000787406 | 4.177316334 | 829515 | 10375 |
| 87 | 4.05216161483e-10 | 7.762981803 | 11.69755744 | 3.93457563699 | 902268 | 10399 |
| 91 | 4.23846789597e-10 | 8.673081643 | 12.512113196 | 3.83903155301 | 943844 | 10385 |
| 92 | 4.28504446625e-10 | 7.685900208 | 11.021693125 | 3.335792917 | 954229 | 10389 |
| 102 | 4.75081016911e-10 | 7.33171763 | 11.297702959 | 3.96598532901 | 1058177 | 10377 |
| 103 | 4.79738673939e-10 | 7.614667959 | 11.020737287 | 3.406069328 | 1068554 | 10386 |


## 9. Instrumented append equivalence

Ten representative real payloads independently matched frozen Journal.append byte-for-byte, with valid Journal.read hash chains, in each of the four contexts. Initial journal bytes and payload family are identical. Setup copying/loading is outside append timing. Secondary atomic JSON also matches exact bytes; NPZ arrays/metadata are validated after each archive.

The starting journal has 1,328,774 bytes; the payload family is ordinary rows
100–109, each approximately 10.4 kB. This is not a toy small-record test.

## 10. Bind burst results

N=120, stalls=0; total append `{"min": 0.004784653996466659, "p50": 0.0057781475043157116, "p75": 0.006835044005129021, "p90": 0.009107438000501135, "p95": 0.012196242045320105, "p99": 0.01354799844964873, "max": 0.014662004992715083}`; non-fsync `{"min": 0.003565562001313083, "p50": 0.004492747000767849, "p75": 0.005172167257114779, "p90": 0.006746118291630413, "p95": 0.01073771409865003, "p99": 0.01223830606948468, "max": 0.012856992994784378}`.

## 11. Bind paced results

60 appends after the first 60 historical ordinary core durations; no delay cap. Sleep total 255.76952902001358 s is outside append statistics. Stall fraction 0.0, fsync statistics `{"min": 0.0013234799989731982, "p50": 0.002252286496513989, "p75": 0.002935340497060679, "p90": 0.004074786597630009, "p95": 0.004687748954893323, "p99": 0.01233082099774032, "max": 0.021658131008734927}`.

## 12. Host results

120 exact burst appends; fsync `{"min": 0.0007924092933535576, "p50": 0.0010048490948975086, "p75": 0.0010893323924392462, "p90": 0.0012944369576871398, "p95": 0.0019359879195690148, "p99": 0.006597637254744774, "max": 0.007344218902289867}`; stalls 0. Only temporary benchmark storage was used.

## 13. Container tmp results

120 exact burst appends; fsync `{"min": 0.0016123800014611334, "p50": 0.0019377790013095364, "p75": 0.002088616747641936, "p90": 0.004963087697979079, "p95": 0.008290234846208477, "p99": 0.008768316450732528, "max": 0.009856613003648818}`; stalls 0. Only temporary benchmark storage was used.

## 14. Docker volume results

120 exact burst appends; fsync `{"min": 0.0015748790028737858, "p50": 0.0019657759985420853, "p75": 0.002139091993740294, "p90": 0.007777740698656999, "p95": 0.008200255447445669, "p99": 0.009393062107264996, "max": 0.013615201998618431}`; stalls 0. Only temporary benchmark storage was used.

## 15. Phase attribution

| phase | p50 [s] | p90 [s] | p95 [s] | max [s] | median stall share |
| --- | --- | --- | --- | --- | --- |
| prepare/hash | 0.000827522511827 | 0.00129292511119 | 0.00158661233436 | 0.00829787401017 | N/A |
| encode | 0.000356262004061 | 0.000631219302886 | 0.000683972200932 | 0.00680413399823 | N/A |
| open | 0.000265728507657 | 0.00041004519444 | 0.0004610401018 | 0.0067068280041 | N/A |
| write | 0.000319740007399 | 0.00056285149476 | 0.000945165799931 | 0.00203924300149 | N/A |
| flush | 1.90200080397e-06 | 2.82299297396e-06 | 3.55779702659e-06 | 1.44169898704e-05 | N/A |
| fsync | 0.00124021049851 | 0.00197484509845 | 0.00295129435253 | 0.00478452300013 | N/A |
| close | 0.00021124050545 | 0.000337434500398 | 0.000398014603707 | 0.00142208499892 | N/A |

All primary bind stall events:

| rep | total [s] | fsync [s] | fsync fraction | file size before | Dirty kB | Writeback kB | CPU [s] | wall [s] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

No events qualify for this table (0/120); median stall-event shares are undefined,
not zero. The original figure retains the 1 s reference. A separately archived
millisecond-scale visual zoom shows the small phase contributions; it changes
neither source measurements nor classification rules.

Frozen checks: `{"five_stalls": false, "median_fsync_share": false, "non_fsync_p95": true, "other_phase_share": false}`. Wilson historical [0.13723724101363852, 0.28673333966375725], bind [0.0, 0.031019166418703486]; overlap=False. Frequency comparison is descriptive only.

## 16. Sync-mode controls

| mode | N | total p50 [s] | p95 [s] | max [s] | stalls |
| --- | --- | --- | --- | --- | --- |
| fsync | 30 | 0.00744723450043 | 0.0177152217926 | 0.019751682994 | 0 |
| fdatasync | 30 | 0.00747501000296 | 0.015910569852 | 0.0178106349922 | 0 |
| flush_only | 30 | 0.00649165950017 | 0.0149961368494 | 0.0254274059989 | 0 |

Temporary controls only; no recommendation to replace production fsync. Optional strace context/status: host; syscall statistics `{"fsync": {"n": 20, "min": 0.000781, "p50": 0.00093, "p75": 0.00104825, "p90": 0.0037994000000000014, "p95": 0.004871300000000001, "p99": 0.006608659999999997, "max": 0.007043}, "fdatasync": {"n": 0, "min": null, "p50": null, "p75": null, "p90": null, "p95": null, "p99": null, "max": null}, "write": {"n": 26, "min": 1.8e-05, "p50": 4.7e-05, "p75": 5.75e-05, "p90": 7.549999999999999e-05, "p95": 8.725e-05, "p99": 0.0006615, "max": 0.000852}, "openat": {"n": 860, "min": 9e-06, "p50": 2.5e-05, "p75": 3.3e-05, "p90": 4.6e-05, "p95": 5.8e-05, "p99": 0.0006422599999999973, "max": 0.006186}, "close": {"n": 789, "min": 8e-06, "p50": 1.6e-05, "p75": 2e-05, "p90": 2.9e-05, "p95": 4.2e-05, "p99": 0.0003834400000000005, "max": 0.006408}, "rename": {"n": 0, "min": null, "p50": null, "p75": null, "p90": null, "p95": null, "p99": null, "max": null}}`.

## 17. Atomic-json secondary results

| context | N | temp-file fsync p50 | p95 | max |
| --- | --- | --- | --- | --- |
| bind | 30 | 0.00230066699442 | 0.00863829315349 | 0.00907848300994 |
| container_tmp | 30 | 0.0087746769932 | 0.0253849727051 | 0.0377010450029 |

Separate distribution; not an ordinary scalar append operation.

## 18. NPZ/archive secondary results

N=10; NPZ fsync `{"min": 0.011912401008885354, "p50": 0.013960023497929797, "p75": 0.01737838000553893, "p90": 0.018488358000467997, "p95": 0.018657202501344727, "p99": 0.01879227810204611, "max": 0.018826047002221458}`; metadata fsync `{"min": 0.0011396310001146048, "p50": 0.0013151810053386725, "p75": 0.0016017895068216603, "p90": 0.0023847735050367192, "p95": 0.0025362052503624, "p99": 0.0026573506466229446, "max": 0.002687636995688081}`; compression/write `{"min": 0.7324231859965948, "p50": 0.8635011760052294, "p75": 0.8828442439989885, "p90": 0.9361055743007455, "p95": 1.0005965066564384, "p99": 1.052189252540993, "max": 1.0650874390121317}`. These costs cannot cause historical F0_C0_A0 latency because that path does not archive NPZ.

## 19. Dirty/writeback correlation

Full stall/nonstall medians, p90 and ranges are in analysis/kernel_state.json. Bind observations: `{"Dirty": {"stall": {"n": 0, "min": null, "p50": null, "p75": null, "p90": null, "p95": null, "p99": null, "max": null}, "nonstall": {"n": 120, "min": 184.0, "p50": 192.0, "p75": 192.0, "p90": 192.0, "p95": 192.0, "p99": 192.0, "max": 192.0}, "units": "kB from proc/meminfo; descriptive, not causal"}, "Writeback": {"stall": {"n": 0, "min": null, "p50": null, "p75": null, "p90": null, "p95": null, "p99": null, "max": null}, "nonstall": {"n": 120, "min": 0.0, "p50": 0.0, "p75": 0.0, "p90": 0.0, "p95": 4.0, "p99": 4.0, "max": 4.0}, "units": "kB from proc/meminfo; descriptive, not causal"}}`. Disk device mapping is used only if directly observed; null means unavailable. Correlation is not causation.

## 20. File-size/cadence analysis

Historical step spacings: `[8, 7, 1, 11, 8, 6, 6, 2, 3, 2, 2, 3, 1, 2, 2, 9, 7, 4, 1, 10, 1]`; bind repetition spacings: `[]`. Correlations: `{"step": -0.06802052749443394, "offset_before": -0.06800489169689258, "line_bytes": -0.03748658751462025, "block": null, "dt": null, "core_s": 0.11170907437043105, "previous_step_core_s": 0.1858338973433463}`. Full offsets, modulo exploratory counts and byte-boundary remainders are retained. No periodicity significance or deterministic size-trigger mechanism is inferred from this sample.

## 21. Operation attribution verdict

JOURNAL FSYNC NOT YET ISOLATED. Phase measurement identifies the present diagnostic operation only; matching historical tail frequency does not retrospectively instrument old syscalls.

In this experiment the historical tail frequency did NOT match. No measured
phase blocked for 7–9 s, so the exact historical blocking phase remains unknown.

## 22. Storage-path attribution verdict

STORAGE PATH UNRESOLVED. Frozen criteria: `{"classification": "STORAGE_PATH_UNRESOLVED", "reason": "operation gate not passed"}`. This is evidence about observed paths, not proof of an OS or Docker defect.

## 23. What is NOT proven

No new PDE result, no full-horizon coupling evidence, no cost authorization, no proof of crash durability under a changed scheme, and no kernel/VM bug established. Host/container controls also differ in execution environment; the path classification is empirical evidence, not elimination of all possible confounders.

Byte-equivalent output does not recreate the producer's original runtime state,
object lifecycle, allocator/GC state, or concurrent host/VM activity. Sleep matches
the elapsed core cadence, not the computational workload. These are limitations,
not evidence that any of those unmeasured mechanisms caused the historical stall.
With no new stalls, Dirty/Writeback stall-versus-nonstall contrasts and benchmark
stall periodicity cannot be estimated.

## 24. Candidate next storage experiments

Recommended next experiment, separately authorized and with no PDE: replay the
two-journal producer sequence (scalar append, then timing append) on host and bind,
using frozen historical payloads and core cadence, unchanged fsync semantics, and
independent host/guest writeback plus Python GC-event observation. This tests whether
the missing producer/storage context matters. It does not assert that the timing
journal caused the old archival latency: its append remains outside archival_s.
Do not migrate storage or reduce durability based on the current result.
Not implemented. Removing fsync is not a harmless optimization.

## 25. Prefix final integrity

Final guard `{"COMPLETE_exists": false, "FAILED_exists": false, "accepted_steps": 127, "administrative_guard_version": "step3a10-numerical-identity-with-separate-session-lineage-v1", "all_inherited_bindings_unchanged": true, "archive_members_identical": 72, "checkpoint_sha256": "bc229cc210b61b6b5ef374f9dd33ed69b88f540065652bc3b415c306e7222d21", "current_session_HEAD": "f78f90ce2c690ab444de5224c3c395ba2f8a4213", "historical_checkpoint_identity_role": "immutable trajectory origin, NOT current session HEAD", "historical_execution_base_HEAD": "5f2d9498e1d0715f57b125ee49b4ddfd533467ae", "interval128_executed": false, "numerical_guard": "PASS", "numerical_identity_sha256": "3c5750c879403b190ec02f9b88592d80348034c649cc48007410ca2061ba090e", "prefix_guard": "PASS", "primary_sessions": 2, "time": 6.148107277665776e-10, "timestamp_UTC": "2026-09-08T15:57:57.724912+00:00"}`. All inherited STEP3A9 prefix and STEP3A10 files remain byte-identical. Accepted=127, no row/field/checkpoint128, no new primary session, no FAILED/COMPLETE. interval128 NOT RUN.

## 26. Tests / CI

Only pure I/O/policy tests; no numerical FEM solve in tests. Results: `{"host_final.xml": {"name": "pytest", "errors": "0", "failures": "0", "skipped": "0", "tests": "57", "time": "0.437", "timestamp": "2026-09-08T18:57:53.798287", "hostname": "victor-ThinkPad-W541"}, "host_preflight.xml": {"name": "pytest", "errors": "0", "failures": "0", "skipped": "0", "tests": "27", "time": "0.264", "timestamp": "2026-09-08T18:28:59.256827", "hostname": "victor-ThinkPad-W541"}, "pinned_preflight.xml": {"name": "pytest", "errors": "0", "failures": "0", "skipped": "0", "tests": "57", "time": "0.554", "timestamp": "2026-09-08T15:31:08.979651+00:00", "hostname": "c0c278bbce6b"}}`. Diagnostic wall 325.48498050961643 s / 3600 s; all batches/cleanup included, pytest/report rendering excluded. Remote CI not observed; no push.

## 27. Exact verdict

JOURNAL FSYNC NOT YET ISOLATED; STORAGE PATH UNRESOLVED. Full CHNS transfer remains incomplete; STEP3A10 authorization=false unchanged. MODEL NOT YET VALIDATED.
