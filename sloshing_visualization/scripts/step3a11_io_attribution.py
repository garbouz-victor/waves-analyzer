"""External storage measurements only. Never constructs a numerical solver."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from datetime import datetime, timezone

import numpy as np
from sloshing.multiphase.phase_rate_history import (
    Journal, atomic_json, archive, load_archive, VERSION,
)
from sloshing.multiphase.linearized_ch import json_hash, array_hash

ROOT = Path("validation_results/step3a11")
PREFIX = Path("validation_results/step3a9/full_coupling_L0")
PRIOR = Path("validation_results/step3a10")
PHASES = ("copy_prepare", "record_hash", "json_encode", "open", "write", "flush", "fsync", "close")


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def output(path, value):
    """New evidence only, never overwrite an experiment."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    atomic_json(path, value)


def utc():
    return datetime.now(timezone.utc).isoformat()


def quantiles(values):
    a = np.asarray(values, dtype=float)
    if not len(a):
        return {k: None for k in ("min", "p50", "p75", "p90", "p95", "p99", "max")}
    if not np.isfinite(a).all():
        raise ValueError("Nonfinite measurements")
    return dict(zip(("min", "p50", "p75", "p90", "p95", "p99", "max"),
                    map(float, np.percentile(a, [0, 50, 75, 90, 95, 99, 100]))))


def wilson(count, n):
    if not n:
        return [None, None]
    z = 1.959963984540054
    p = count/n; den = 1+z*z/n
    mid = (p+z*z/(2*n))/den
    radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0., mid-radius), min(1., mid+radius)]


def command(args):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=20)
        return {"argv": args, "returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": args, "not_available": str(exc)}


def environment(path):
    path = Path(path).resolve()
    raw = Path("/proc/self/mountinfo").read_text() if Path("/proc/self/mountinfo").exists() else ""
    matches = []
    for line in raw.splitlines():
        left, right = line.split(" - ", 1); a = left.split(); b = right.split()
        mount = a[4].replace("\\040", " ")
        if path == Path(mount) or Path(mount) in path.parents:
            matches.append({"mountpoint": mount, "device": a[2], "options": a[5],
                            "filesystem": b[0], "source": b[1], "super_options": b[2]})
    mount = max(matches, key=lambda m: len(m["mountpoint"])) if matches else None
    st = path.stat()
    return {"path": str(path), "uid": os.getuid(), "gid": os.getgid(), "uname": list(os.uname()),
            "python": sys.version, "device_id": st.st_dev, "device_major_minor": [os.major(st.st_dev), os.minor(st.st_dev)],
            "mount": mount, "findmnt": command(["findmnt", "-J", "-T", str(path)]),
            "stat_f": command(["stat", "-f", str(path)]), "mountinfo": raw,
            "strace": shutil.which("strace")}


def proc_state(device=None):
    result = {}
    selections = {"meminfo": ("Dirty", "Writeback", "WritebackTmp", "Buffers", "Cached"),
                  "vmstat": ("nr_dirty", "nr_writeback"),
                  "self/io": ("wchar", "write_bytes", "cancelled_write_bytes")}
    for name, keys in selections.items():
        vals = {k: None for k in keys}
        try:
            for line in Path("/proc", name).read_text().splitlines():
                parts = line.replace(":", "").split()
                if parts and parts[0] in vals:
                    vals[parts[0]] = int(parts[1])
        except (OSError, ValueError):
            pass
        result[name] = vals
    result["diskstats"] = None
    if device is not None:
        try:
            for line in Path("/proc/diskstats").read_text().splitlines():
                p = line.split()
                if (int(p[0]), int(p[1])) == tuple(device):
                    result["diskstats"] = {"device": p[2], "counters": list(map(int, p[3:]))}
        except (OSError, ValueError):
            pass
    return result


class Timer:
    def __init__(self):
        self.wall = time.perf_counter(); self.cpu = time.process_time(); self.values = {}

    def call(self, name, fn):
        w = time.perf_counter(); c = time.process_time()
        value = fn()
        self.values[name+"_s"] = time.perf_counter()-w
        self.values[name+"_cpu_s"] = time.process_time()-c
        return value

    def finish(self, name):
        self.values[name+"_s"] = time.perf_counter()-self.wall
        self.values[name+"_cpu_s"] = time.process_time()-self.cpu
        return self.values


def sync_call(timer, stream, mode, probes, device):
    before = timer.call("probe_before", lambda: proc_state(device)) if probes else None
    if mode == "fsync":
        timer.call("fsync", lambda: os.fsync(stream.fileno()))
    elif mode == "fdatasync":
        timer.call("fsync", lambda: os.fdatasync(stream.fileno()))
    elif mode == "flush_only":
        timer.values.update(fsync_s=0., fsync_cpu_s=0.)
    else:
        raise ValueError(mode)
    after = timer.call("probe_after", lambda: proc_state(device)) if probes else None
    timer.values["kernel_before_fsync"] = before
    timer.values["kernel_after_fsync"] = after
    if before and before["diskstats"] and after and after["diskstats"]:
        timer.values["diskstats_delta"] = [b-a for a,b in zip(before["diskstats"]["counters"], after["diskstats"]["counters"])]
    else:
        timer.values["diskstats_delta"] = None


class InstrumentedJournal:
    def __init__(self, path):
        self.path = Path(path)
        self.rows = Journal.read(self.path) if self.path.exists() else []

    def append(self, value, *, mode="fsync", probes=True, repetition=0):
        size_before = self.path.stat().st_size if self.path.exists() else 0
        dev = self.path.parent.stat().st_dev
        timer = Timer(); stamp = timer.wall
        row = timer.call("copy_prepare", lambda: {**copy.deepcopy(value),
            "previous_record_sha256": self.rows[-1]["record_sha256"] if self.rows else None})
        row["record_sha256"] = timer.call("record_hash", lambda: json_hash(row))
        stream = timer.call("open", lambda: self.path.open("ab"))
        try:
            data = timer.call("json_encode", lambda: (json.dumps(row, sort_keys=True, allow_nan=False)+"\n").encode())
            timer.call("write", lambda: stream.write(data))
            timer.call("flush", stream.flush)
            sync_call(timer, stream, mode, probes, (os.major(dev), os.minor(dev)))
        finally:
            timer.call("close", stream.close)
        timer.call("rows_publish", lambda: self.rows.append(row))
        rec = timer.finish("total_append")
        rec.update(repetition=repetition, monotonic_timestamp=stamp, mode=mode,
                   line_bytes=len(data), file_size_before=size_before,
                   file_size_after=self.path.stat().st_size)
        rec["total_io_s"] = sum(rec[k+"_s"] for k in ("open", "write", "flush", "fsync", "close"))
        rec["encode_hash_s"] = sum(rec[k+"_s"] for k in ("copy_prepare", "record_hash", "json_encode"))
        rec["non_fsync_s"] = rec["total_append_s"]-rec["fsync_s"]
        rec["stall"] = rec["total_append_s"] >= 1.0
        rec["phase_fractions"] = {k: rec[k+"_s"]/rec["total_append_s"] for k in PHASES}
        rec["encode_fraction"] = rec["encode_hash_s"]/rec["total_append_s"]
        return row, rec


def instrumented_atomic(path, value, *, probes=True):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    timer = Timer()
    data = timer.call("json_encode", lambda: (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode())
    fd, name = timer.call("mkstemp", lambda: tempfile.mkstemp(prefix=".pending-", dir=path.parent))
    stream = timer.call("fdopen", lambda: os.fdopen(fd, "wb"))
    try:
        timer.call("write", lambda: stream.write(data)); timer.call("flush", stream.flush)
        sync_call(timer, stream, "fsync", probes, None)
    finally:
        timer.call("close", stream.close)
    timer.call("replace", lambda: os.replace(name, path))
    result = timer.finish("total_atomic"); result["file_bytes"] = len(data)
    return result


def instrumented_archive(path, arrays, metadata):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists(): raise FileExistsError(path)
    timer = Timer()
    arrays = timer.call("array_copy", lambda: {k: np.array(v, copy=True) for k,v in arrays.items()})
    if metadata.get("accepted") is not False and not all(np.isfinite(v).all() for v in arrays.values()):
        raise ValueError("Nonfinite archive")
    record = timer.call("array_hash", lambda: {**metadata, "archive_version": VERSION,
        "array_sha256": {k: array_hash(v) for k,v in arrays.items()}})
    record["semantic_sha256"] = timer.call("metadata_hash", lambda: json_hash(record))
    pending = Path(timer.call("mkdtemp", lambda: tempfile.mkdtemp(prefix=".pending-", dir=path.parent)))
    stream = timer.call("open", lambda: (pending/"arrays.npz").open("wb"))
    try:
        timer.call("npz_compress_write", lambda: np.savez_compressed(stream, **arrays))
        timer.call("flush", stream.flush); sync_call(timer, stream, "fsync", True, None)
    finally:
        timer.call("close", stream.close)
    meta_measure = timer.call("metadata_atomic", lambda: instrumented_atomic(pending/"metadata.json", record))
    timer.call("rename", lambda: os.rename(pending, path))
    result = timer.finish("total_archive")
    result["metadata_atomic"] = meta_measure
    result["npz_bytes"] = (path/"arrays.npz").stat().st_size
    return result


def inherited_files():
    paths = sorted(p for folder in (PREFIX, PRIOR) for p in folder.rglob("*") if p.is_file())
    return {str(p): {"sha256": sha(p), "bytes": p.stat().st_size} for p in paths}


def check_bindings(bindings):
    for p, info in bindings.items():
        if sha(p) != info["sha256"] or Path(p).stat().st_size != info["bytes"]:
            raise ValueError("Inherited file altered: "+p)


def guard(final=False):
    from step3a10_admin import numerical_guard, prefix_guard, lineage
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    f = numerical_guard(); manifest = prefix_guard(); members = 0
    with tarfile.open(d.ROOT/"policy/production_source.tar.gz") as tar:
        for member in tar.getmembers():
            if member.isfile():
                p = d.SOURCE/member.name[len("multiphase/"):] if member.name.startswith("multiphase/") else Path(member.name)
                if p.read_bytes() != tar.extractfile(member).read(): raise ValueError("Archive member mismatch")
                members += 1
    if manifest["checkpoint_semantic_sha256"] != "bc229cc210b61b6b5ef374f9dd33ed69b88f540065652bc3b415c306e7222d21":
        raise ValueError("Wrong checkpoint127")
    assert len(list(PREFIX.glob("sessions/*/session_status.json"))) == 2
    assert not (PREFIX/"fields/000128").exists() and not (PREFIX/"checkpoints/000128").exists()
    bindings = inherited_files()
    if final:
        check_bindings(read(ROOT/"preflight/prefix_integrity.json")["files"])
        if set(bindings) != set(read(ROOT/"preflight/prefix_integrity.json")["files"]):
            raise ValueError("New inherited files/session detected")
    result = {"timestamp_UTC": utc(), **lineage(f), "numerical_guard": "PASS", "prefix_guard": "PASS",
        "archive_members_identical": members, "accepted_steps": 127, "interval128_executed": False,
        "time": manifest["time"], "checkpoint_sha256": manifest["checkpoint_semantic_sha256"],
        "primary_sessions": 2, "FAILED_exists": False, "COMPLETE_exists": False}
    if final:
        result["all_inherited_bindings_unchanged"] = True
        output(ROOT/"final_integrity.json", result)
    else:
        output(ROOT/"preflight/numerical_guard.json", result)
        output(ROOT/"preflight/prefix_integrity.json", {**result, "files": bindings})
    print(json.dumps(result, indent=2))


def historical():
    check_bindings(read(ROOT/"preflight/prefix_integrity.json")["files"])
    rows = Journal.read(PREFIX/"scalar_history.jsonl")
    timings = {r["step"]: r["timing"] for r in Journal.read(PREFIX/"timing_history.jsonl")}
    events = read(PRIOR/"cost_audit/event_schedule.json")["events"]
    offsets = []; offset = 0
    for line in (PREFIX/"scalar_history.jsonl").read_bytes().splitlines(keepends=True):
        row = json.loads(line); offsets.append({"step": row["step"], "offset_before": offset,
            "line_bytes": len(line), "file_size_after": offset+len(line)})
        offset += len(line)
    ordinary = []
    for e in events:
        i = e["step"]
        if i > 127 or e["archive_class"] != "F0_C0_A0": continue
        t = timings[i]; core = t["total_s"]-t["archival_s"]
        if core < 0: raise ValueError("Invalid historical core")
        ordinary.append({**e, **offsets[i], "total_s": t["total_s"], "core_s": core,
            "archival_s": t["archival_s"], "stall": t["archival_s"] >= 1.,
            "previous_step_core_s": timings[i-1]["total_s"]-timings[i-1]["archival_s"] if i > 1 else None})
    stalls = [r for r in ordinary if r["stall"]]
    result = {"n": len(ordinary), "statistics": quantiles([r["archival_s"] for r in ordinary]),
        "stall_count": len(stalls), "stall_fraction": len(stalls)/len(ordinary), "wilson95": wilson(len(stalls), len(ordinary)),
        "stall_steps": [r["step"] for r in stalls], "stall_step_spacing": np.diff([r["step"] for r in stalls]).tolist(),
        "records": ordinary, "historical_wall_timestamps": "not recorded per append; physical time is not wall time"}
    output(ROOT/"historical/ordinary_archive.json", result)
    output(ROOT/"historical/journal_offsets.json", offsets)
    return rows, ordinary


def policy_guard():
    p = read(ROOT/"policy.json")
    if sha(ROOT/"policy.json") != (ROOT/"policy_sha256").read_text().strip(): raise ValueError("Policy changed")
    for name, digest in p["source_sha256"].items():
        if sha(name) != digest: raise ValueError("Attribution implementation changed: "+name)
    return p


def payloads(policy):
    rows = Journal.read(PREFIX/"scalar_history.jsonl")
    return [{k:v for k,v in rows[i].items() if k not in ("record_sha256", "previous_record_sha256")}
            for i in policy["payload_steps"]]


def safe_cleanup(path, parent):
    path = Path(path).resolve(); parent = Path(parent).resolve()
    if path.parent != parent or not path.name.startswith("step3a11-") or not (path/"STEP3A11_TEMP.json").exists():
        raise ValueError("Unsafe temporary cleanup")
    shutil.rmtree(path)


def equivalence(folder, family, initial):
    a = folder/"equivalence-production.jsonl"; b = folder/"equivalence-instrumented.jsonl"
    a.write_bytes(initial); b.write_bytes(initial)  # setup, outside append timers
    production = Journal(a); replica = InstrumentedJournal(b); results = []
    for i, value in enumerate(family, 1):
        production.append(value); replica.append(value, repetition=i)
        if a.read_bytes() != b.read_bytes(): raise ValueError("Replica is not byte-equivalent")
        assert Journal.read(a) == Journal.read(b)
        results.append({"payload_index": i, "sha256": sha(a), "bytes": a.stat().st_size, "byte_equal": True})
    return {"passed": True, "n": len(results), "records": results}


def operation_classification(records):
    stalls = [r for r in records if r["total_append_s"] >= 1.]
    med = lambda key: float(np.median([r[key]/r["total_append_s"] for r in stalls])) if stalls else None
    fs = med("fsync_s")
    other = {k: med(k+"_s") for k in PHASES if k != "fsync"}
    checks = {"five_stalls": len(stalls) >= 5, "median_fsync_share": fs is not None and fs >= .90,
              "non_fsync_p95": bool(records) and quantiles([r["non_fsync_s"] for r in records])["p95"] <= .10,
              "other_phase_share": bool(stalls) and all(v < .20 for v in other.values())}
    return {"classification": "FSYNC_DOMINATED" if all(checks.values()) else "FSYNC_NOT_YET_ISOLATED",
            "checks": checks, "stall_count": len(stalls), "median_fsync_share": fs, "other_median_shares": other}


def context_classification(metrics, operation_pass):
    if not operation_pass: return {"classification": "STORAGE_PATH_UNRESOLVED", "reason": "operation gate not passed"}
    a,b,c,d = [metrics[k] for k in ("bind", "host", "container_tmp", "docker_volume")]
    ps = [v["p95"] for v in (a,b,c,d)]; fs = [v["stall_fraction"] for v in (a,b,c,d)]
    bind = ps[0] >= 5*max(*ps[1:],1e-6) and fs[0] >= .05 and sum(f <= .02 for f in fs[1:]) >= 2
    host = (min(fs[:2]) >= .05 and max(ps[:2])/max(min(ps[:2]),1e-30) < 2 and
            any(v["stall_fraction"] <= .02 and v["p95"]*5 <= min(ps[:2]) for v in (c,d)))
    docker = (min(fs[i] for i in (0,2,3)) >= .05 and fs[1] <= .02 and
              all(ps[i] >= 5*max(ps[1],1e-6) for i in (0,2,3)))
    return {"classification": "BIND_MOUNT_SPECIFIC_EVIDENCE" if bind else "HOST_FILESYSTEM_DURABILITY_EVIDENCE" if host else
            "DOCKER_STORAGE_PATH_EVIDENCE" if docker else "STORAGE_PATH_UNRESOLVED",
            "checks": {"bind_specific": bind, "host_filesystem": host, "docker_storage": docker},
            "metrics": metrics, "bind_over_controls_p95": ps[0]/max(*ps[1:],1e-6)}


def worker(kind, context, base, budget_s):
    p = policy_guard(); family = payloads(p); started = time.monotonic(); deadline = started+budget_s
    base = Path(base).resolve()
    temp = Path(tempfile.mkdtemp(prefix="step3a11-", dir=base))
    (temp/"STEP3A11_TEMP.json").write_text(json.dumps({"owner": "step3a11", "kind": kind}))
    result = {"kind": kind, "context": context, "start_UTC": utc(), "environment": environment(temp),
              "policy_sha256": sha(ROOT/"policy.json"), "records": [], "status": "running"}
    def check(delay=0.):
        if time.monotonic()+delay >= deadline: raise TimeoutError("Frozen diagnostic budget exhausted")
    try:
        initial = (PREFIX/"scalar_history.jsonl").read_bytes()
        result["starting_journal_sha256"] = hashlib.sha256(initial).hexdigest()
        result["starting_journal_bytes"] = len(initial)
        if kind == "burst":
            check(); result["equivalence"] = equivalence(temp, family, initial)
        if kind in ("burst", "paced", "fsync", "fdatasync", "flush_only", "trace"):
            if kind == "fdatasync" and not hasattr(os, "fdatasync"):
                result["status"] = "not_available"
            else:
                n = p["counts"][kind] if kind in ("burst", "paced", "trace") else p["counts"]["sync_control"]
                path = temp/"journal.jsonl"; path.write_bytes(initial)
                journal = InstrumentedJournal(path)
                for i in range(n):
                    delay = p["paced_delays"][i]["core_s"] if kind == "paced" else 0.
                    check(delay)
                    if delay: time.sleep(delay)
                    _, rec = journal.append(family[i%len(family)], mode=kind if kind in ("fdatasync", "flush_only") else "fsync", repetition=i+1)
                    rec["pacing_sleep_s"] = delay; result["records"].append(rec)
                Journal.read(path)
                result.update(expected_n=n, status="complete", final_journal_sha256=sha(path), final_journal_bytes=path.stat().st_size)
        elif kind == "atomic":
            # Verify the secondary replica too; setup and equivalence excluded from N.
            atomic_json(temp/"original.json", family[0]); instrumented_atomic(temp/"replica.json", family[0])
            assert (temp/"original.json").read_bytes() == (temp/"replica.json").read_bytes()
            for i in range(p["counts"]["atomic"]):
                check(); rec = instrumented_atomic(temp/"atomic.json", family[i%len(family)])
                result["records"].append({"repetition": i+1, **rec})
            result.update(status="complete", expected_n=p["counts"]["atomic"], byte_equivalence=True)
        elif kind == "archive":
            arrays, meta = load_archive(PREFIX/"checkpoints/000127")
            meta = {k:v for k,v in meta.items() if k not in ("array_sha256", "semantic_sha256", "archive_version")}
            for i in range(p["counts"]["archive"]):
                check(); target = temp/f"checkpoint-copy-{i:02d}"
                rec = instrumented_archive(target, arrays, meta)
                loaded, metadata = load_archive(target)
                assert all(np.array_equal(loaded[k], v) for k,v in arrays.items())
                rec.update(repetition=i+1, semantic_sha256=metadata["semantic_sha256"], npz_sha256=sha(target/"arrays.npz"))
                result["records"].append(rec)
            result.update(status="complete", expected_n=p["counts"]["archive"])
        else: raise ValueError(kind)
        check()
    except Exception as exc:
        result.update(status="incomplete", error=repr(exc))
    finally:
        # Export detailed measurements before removing any temporary benchmark file.
        result["payload_description_sha256"] = json_hash({"initial": result.get("starting_journal_sha256"), "payload_steps": p["payload_steps"]})
        result["wall_before_cleanup_s"] = time.monotonic()-started
        print(json.dumps({"measurements_before_cleanup": result}, allow_nan=False), flush=True)
        safe_cleanup(temp, base)
        result.update(temporary_files_removed=True, wall_s=time.monotonic()-started)
    print(json.dumps({"final_result": result}, allow_nan=False), flush=True)
    return result["status"] in ("complete", "not_available")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("preflight", "final", "historical", "worker"))
    parser.add_argument("--kind"); parser.add_argument("--context"); parser.add_argument("--base", default=".")
    parser.add_argument("--budget-s", type=float, default=3585.)
    args = parser.parse_args()
    if args.action in ("preflight", "final"): guard(args.action == "final")
    elif args.action == "historical": historical()
    else: raise SystemExit(0 if worker(args.kind, args.context, args.base, args.budget_s) else 2)
