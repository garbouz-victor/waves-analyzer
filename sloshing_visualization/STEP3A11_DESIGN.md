# STEP3A11 — journal/fsync latency attribution (no PDE)

Execution base: f78f90ce2c690ab444de5224c3c395ba2f8a4213.
The scientific prefix remains 127/1068, t=6.148107277665776e-10 s.
STEP3A10 authorization=false and its 14400/19800 s caps are immutable history.
No continuation authorization, numerical solve, storage optimization, or source
change under src/sloshing/multiphase is permitted in this iteration.

## Source audit

RateTrajectory.advance starts its archival timer immediately before scalar
Journal.append(row). Selected residual JSON, field, and checkpoint operations
follow conditionally. Both archival_s and total_s stop BEFORE timing_journal.append.
FullRateTrajectory disables the isolated expanded-audit file path; full constraint
audits are in observe(), before archival. DivergenceTrajectory inherits advance.
For F0_C0_A0, only scalar Journal.append is inside archival_s. Neither NPZ nor
timing_history append can explain a stall in this particular measured interval.

Exact append order is deepcopy/previous hash, json_hash, open(ab), JSON serialization
and encode, buffered write, flush, fsync, close, in-memory rows.append. In particular,
JSON encoding occurs AFTER open in the actual frozen implementation. The external
replica preserves that order, formatting, hash computation and file buffering.
File copying, directory creation, journal loading, and file-size stat calls are
outside append timers. All per-append results stay in memory until the batch ends.

Atomic_json fsyncs its temporary file before os.replace. Archive compresses arrays,
fsyncs the NPZ, uses atomic_json for metadata, then renames the directory. Neither
implementation fsyncs the parent directory. This is a source observation, not a
change: atomic publication and crash durability are distinct guarantees.

## Frozen experiment

The machine policy is written and hashed before scientific timing benchmarks.
Preflight reuses numerical_guard and prefix_guard without modifying their code or
environment identity. Archive members and all inherited prefix/STEP3A10 file bytes
are checked before and after. No solver object is constructed.

All exact append contexts start with a byte copy of the complete scalar prefix
0..127. Payloads cycle through the first ten ordinary rows with step>=100, removing
both hash fields before natural append reconstruction. Ten representative appends
are checked against frozen Journal.append for exact byte equality and valid chains
in each context before using its instrumentation. Equivalence is not included in
the N measured primary appends. A fresh copy is used for every batch.

Fixed sequential order: bind burst (120), host burst (120), container /tmp burst
(120), unique Docker volume burst (120), bind paced (60), bind fsync/fdatasync/
flush-only controls (30 each), atomic_json bind/tmp (30 each), checkpoint127 NPZ
archive on bind (10), optional strace sequence (20). No discarded or adaptive
repetitions. Paced delays are the first 60 historical ordinary core_s observations
in chronological order, with no cap; their exact steps/delays are frozen in policy.
No sleep duration enters append latency. Docker context/image remain production's;
/work is the same repository bind. Other contexts use only temporary files, not a
moved scientific trajectory. Strace is used only if already installed, preferring
container bind, otherwise host; absence is reported without package installation.
An installed but unusable optional strace is recorded as unavailable after one
attempt, without retrying in a different context or invalidating mandatory probes.

## Instrumentation and rules

Per phase: wall and process CPU time. Immediately around fsync, read available
Dirty/Writeback/WritebackTmp/Buffers/Cached, nr_dirty/nr_writeback and process I/O.
Probe wall time is reported separately AND included in total_append and conservative
non_fsync_s=total_append_s-fsync_s. fsync_s brackets only the durability call.
Disk counters are used only when stat device major/minor directly matches diskstats;
otherwise not_available, with no Docker/VM device inference. No global sync/drop-cache.

Stall means total append (or historical archival) >=1.0 s, descriptively only.
Primary operation decision uses ONLY bind burst's 120 exact fsync appends, never
pooling paced/control samples to reach the five-stall minimum. FSYNC_DOMINATED needs
>=5 stalls, median fsync share>=0.90, p95(total minus fsync)<=0.10 s, and every other
single measured phase's median stall share<0.20. Otherwise FSYNC_NOT_YET_ISOLATED.

Path classification is attempted only after operation qualification. Bind-specific:
bind fsync p95>=5*max(host,tmp,volume p95,1e-6), bind stall fraction>=0.05 and at least
two controls with fraction<=0.02. Host-filesystem: host and bind fractions>=0.05,
their fsync p95 ratio<2, and tmp or volume has fraction<=0.02 and fsync p95 at least
5x smaller than min(host,bind). Docker-storage: bind/tmp/volume fractions>=0.05,
host fraction<=0.02, and each Docker-context fsync p95>=5*max(host p95,1e-6).
These explicit 5x/0.02 definitions resolve 'materially faster' before measurement.
No classification match means unresolved; no claim of an OS/Docker bug.

Percentiles use deterministic linear interpolation. Wilson 95% intervals use
z=1.959963984540054. Frequency overlap, size/cadence correlations, mod 2/4/5/8/10/16,
and Dirty/Writeback comparisons are diagnostic only. Byte-boundary remainders use
4 KiB, 64 KiB and 1 MiB descriptively. No outlier deletion or fitted
classification threshold. Historical monotonic timestamps are not invented.

Total diagnostic wall cap is 3600 s from start of the sequential benchmark suite,
including equivalence, setup, pacing, process startup, result export, and cleanup.
Normal pytest, policy/preflight inspection and report rendering are excluded.
Reserve 15 s for orderly cleanup; exhaustion produces incomplete batches and STOP,
not a shortened successful N or changed budget. A blocking syscall already entered
cannot be promised interruptible; any overrun is recorded explicitly.

## Outcomes

Even a positive operation/path attribution is an execution/storage finding only.
No cost reauthorization, no interval128, no reduction of fsync, no storage migration.
Future durability changes require separate byte-equivalence and crash/restart
qualification. Overall MODEL NOT YET VALIDATED; full CHNS transfer incomplete.
