"""Administrative host controller for temporary I/O experiments; no PDE path."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from step3a11_io_attribution import (ROOT, PRIOR, PREFIX, output, read, sha, utc,
    historical, environment, policy_guard, check_bindings)


def run_command(args, **kwargs):
    return subprocess.run(args, check=True, capture_output=True, text=True, **kwargs).stdout.strip()


def freeze():
    if (ROOT/"policy.json").exists(): raise FileExistsError("Policy already frozen")
    assert read(ROOT/"preflight/numerical_guard.json")["numerical_guard"] == "PASS"
    rows, ordinary = historical()
    context = run_command(["docker", "context", "show"])
    digest = run_command(["docker", "image", "inspect", "--format", "{{.Id}}", "waves-step3:0.10.0"])
    info = json.loads(run_command(["docker", "info", "--format", "{{json .}}"] ))
    version = json.loads(run_command(["docker", "version", "--format", "{{json .}}"] ))
    host = environment(Path.cwd())
    host.update(docker_context=context, image_digest=digest, docker_version=version,
        docker_storage={k:info.get(k) for k in ("Driver", "DriverStatus", "DockerRootDir", "KernelVersion",
            "OperatingSystem", "ServerVersion", "NCPU", "MemTotal", "ContainersRunning")})
    output(ROOT/"environment/host.json", host)
    names = ["STEP3A11_DESIGN.md", "scripts/step3a11_io_attribution.py", "scripts/step3a11_host_io.py",
             "scripts/step3a11_report.py", "tests/test_step3a11_io.py"]
    policy = {"version": "step3a11-latency-attribution-v1", "timestamp_UTC": utc(),
        "base_HEAD": run_command(["git", "rev-parse", "HEAD"]), "docker_context": context, "image_digest": digest,
        "source_sha256": {p:sha(p) for p in names}, "stall_s": 1.0,
        "counts": {"burst": 120, "paced": 60, "sync_control": 30, "atomic": 30, "archive": 10, "trace": 20},
        "context_order": ["bind", "host", "container_tmp", "docker_volume"],
        "paced_delays": [{"step":r["step"], "core_s":r["core_s"]} for r in ordinary[:60]],
        "pacing_cap_s": None, "payload_steps": [r["step"] for r in ordinary if r["step"] >= 100][:10],
        "starting_journal_sha256": sha(PREFIX/"scalar_history.jsonl"),
        "starting_journal_bytes": (PREFIX/"scalar_history.jsonl").stat().st_size,
        "operation": {"primary": "bind_burst ONLY", "min_stalls": 5, "median_fsync_fraction_min": .90,
            "non_fsync_p95_max_s": .10, "other_phase_median_fraction_strict_max": .20},
        "path": {"bind_p95_ratio_min": 5., "slow_stall_fraction_min": .05, "fast_stall_fraction_max": .02,
            "required_fast_controls": 2, "host_bind_ratio_strict_max": 2., "materially_faster_ratio": 5., "p95_floor_s": 1e-6},
        "budget_s": 3600., "cleanup_reserve_s": 15., "phase_timer": "perf_counter wall + process_time CPU",
        "kernel_probes": "before/after durability call; overhead included in total/non_fsync, not fsync",
        "percentile": "linear interpolation", "wilson_z": 1.959963984540054,
        "exploratory_byte_boundaries": [4096,65536,1048576],
        "no_PDE": True, "no_continuation_authorization": True, "no_production_storage_changes": True}
    assert len(policy["payload_steps"]) == 10 and len(policy["paced_delays"]) == 60
    output(ROOT/"policy.json", policy)
    (ROOT/"policy_sha256").write_text(sha(ROOT/"policy.json")+"\n")
    print(json.dumps({"policy_sha256":sha(ROOT/"policy.json"), "payload_steps":policy["payload_steps"],
                      "pacing_total_s":sum(r["core_s"] for r in policy["paced_delays"]), "docker_context":context}, indent=2))


def parse_worker(stdout):
    records = [json.loads(line) for line in stdout.splitlines() if line.startswith("{")]
    matches = [r["final_result"] for r in records if "final_result" in r]
    if len(matches) != 1: raise ValueError("Worker did not export one final result")
    return matches[0]


def run_suite():
    p = policy_guard()
    check_bindings(read(ROOT/"preflight/prefix_integrity.json")["files"])
    if (ROOT/"benchmarks/suite_status.json").exists(): raise FileExistsError("No automatic experiment rerun")
    volume = "step3a11-"+uuid.uuid4().hex
    project = Path.cwd().resolve(); started = time.monotonic(); deadline = started+p["budget_s"]
    env = os.environ.copy(); env["PYTHONPATH"] = str(project/"src")+os.pathsep+env.get("PYTHONPATH", "")
    created = False; completed = []; status = {"status":"running", "start_UTC":utc(), "volume":volume, "policy_sha256":sha(ROOT/"policy.json")}
    def remaining():
        n = deadline-time.monotonic()-p["cleanup_reserve_s"]
        if n <= 0: raise TimeoutError("Frozen diagnostic wall budget exhausted")
        return n
    def docker_base(context, name):
        rw = context == "bind"
        args = ["docker", "run", "--rm", "--network=none", "--name", name,
            "--user", "0:0" if p["docker_context"] == "desktop-linux" else f"{os.getuid()}:{os.getgid()}",
            "--mount", f"type=bind,source={project},target=/work"+("" if rw else ",readonly"),
            "--env", "XDG_CACHE_HOME=/tmp/step3a11-cache", "--env", "MPLCONFIGDIR=/tmp/step3a11-matplotlib"]
        if context == "docker_volume": args += ["--mount", f"type=volume,source={volume},target=/bench"]
        return args+[p["image_digest"]]
    def one(kind, context, filename, trace=False):
        if run_command(["docker", "context", "show"]) != p["docker_context"]: raise ValueError("Docker context changed")
        name = "step3a11-worker-"+uuid.uuid4().hex[:12]
        base = "." if context in ("bind", "host") else "/tmp" if context == "container_tmp" else "/bench"
        py = ["python3" if context != "host" else sys.executable, "scripts/step3a11_io_attribution.py", "worker",
              "--kind", kind, "--context", context, "--base", base, "--budget-s", str(remaining())]
        trace_path = ROOT/"benchmarks/strace.log"
        if trace:
            py = ["strace", "-ttT", "-e", "trace=openat,write,fsync,fdatasync,close,rename", "-o", str(trace_path)]+py
        args = py if context == "host" else docker_base(context, name)+py
        before = time.monotonic()
        try:
            proc = subprocess.run(args, capture_output=True, text=True, timeout=remaining(), env=env if context == "host" else None)
        except subprocess.TimeoutExpired:
            if context != "host": run_command(["docker", "stop", "-t", "5", name])
            raise
        result = parse_worker(proc.stdout)
        result.update(host_observed_wall_s=time.monotonic()-before, command=args, returncode=proc.returncode, stderr=proc.stderr)
        output(ROOT/"benchmarks"/filename, result)
        completed.append(filename)
        if kind == "burst":
            if not (ROOT/"environment"/(context+".json")).exists(): output(ROOT/"environment"/(context+".json"), result["environment"])
            output(ROOT/"equivalence"/(context+".json"), result["equivalence"])
        print(f"{context}/{kind}: {result['status']}, N={len(result['records'])}, wall={result['host_observed_wall_s']:.3f}s", flush=True)
        if proc.returncode or result["status"] not in ("complete", "not_available"):
            raise RuntimeError("Incomplete benchmark; no retries: "+filename)
        return result
    try:
        # All filesystem writes concern uniquely named temporary benchmark storage.
        actual = run_command(["docker", "volume", "create", "--label", "owner=step3a11", "--label", "purpose=temporary-io-attribution", volume])
        assert actual == volume; created = True
        output(ROOT/"environment/volume_creation.json", json.loads(run_command(["docker", "volume", "inspect", volume])))
        bursts = {}
        for context in p["context_order"]:
            bursts[context] = one("burst", context, context+"_burst.json")
        one("paced", "bind", "bind_paced.json")
        for mode in ("fsync", "fdatasync", "flush_only"): one(mode, "bind", "control_"+mode+".json")
        one("atomic", "bind", "atomic_bind.json"); one("atomic", "container_tmp", "atomic_container_tmp.json")
        one("archive", "bind", "archive_npz.json")
        trace_context = "bind" if bursts["bind"]["environment"]["strace"] else "host" if shutil.which("strace") else None
        if trace_context:
            try:
                one("trace", trace_context, "strace_sequence.json", trace=True)
            except (ValueError, RuntimeError) as exc:
                # Optional diagnostic may lack ptrace permission; no retry elsewhere.
                if not (ROOT/"benchmarks/strace_sequence.json").exists():
                    output(ROOT/"benchmarks/strace_sequence.json", {"status":"NOT AVAILABLE", "error":repr(exc), "context":trace_context})
                status["optional_strace_error"] = repr(exc)
        else: output(ROOT/"benchmarks/strace_sequence.json", {"status":"NOT AVAILABLE"})
        status["status"] = "complete"
    except Exception as exc:
        status.update(status="incomplete", error=repr(exc))
    finally:
        if created:
            inspection = json.loads(run_command(["docker", "volume", "inspect", volume]))[0]
            assert volume.startswith("step3a11-") and inspection["Labels"].get("owner") == "step3a11"
            # Only delete our unique volume after results are exported and it is empty.
            empty = run_command(["docker", "run", "--rm", "--network=none", "--mount",
                f"type=volume,source={volume},target=/bench,readonly", p["image_digest"], "python3", "-c",
                "import os,json; print(json.dumps(os.listdir('/bench')))"])
            output(ROOT/"environment/volume_before_cleanup.json", {"inspection":inspection, "remaining_files":json.loads(empty), "results_exported":completed})
            if json.loads(empty):
                status.update(status="incomplete", cleanup="volume retained: nonempty")
            else:
                run_command(["docker", "volume", "rm", volume]); status["volume_removed_after_empty_check"] = True
        status.update(wall_s=time.monotonic()-started, completed_batches=completed,
                      cap_s=p["budget_s"], actual_cap_pass=time.monotonic()-started<=p["budget_s"], end_UTC=utc())
        output(ROOT/"benchmarks/suite_status.json", status)
    print(json.dumps(status, indent=2))
    return status["status"] == "complete" and status["actual_cap_pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("freeze", "run"))
    args = parser.parse_args()
    if args.action == "freeze": freeze()
    else: raise SystemExit(0 if run_suite() else 2)
