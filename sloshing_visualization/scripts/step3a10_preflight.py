"""Read-only STEP3A10 first hard gate; no cost model or PDE execution.

The original source_guard is called without replacing it or its environment.
A failed guard allows only diagnostic hash checks and a new stop checkpoint.
Writes are confined to new STEP3A10 artifacts, never the inherited trajectory.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import tarfile

ROOT = Path("validation_results/step3a10")
BASE = "c718347d81925faa9a6a8184f4dd4c2d0298aaae"
EXPECTED_SCHEDULE = "21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bindings(paths):
    return {str(p): {"sha256": digest(p), "bytes": p.stat().st_size}
            for p in sorted(map(Path, paths))}


def verify_bindings(records):
    for path, record in records.items():
        if digest(path) != record["sha256"] or Path(path).stat().st_size != record["bytes"]:
            raise ValueError("Inherited file changed: " + path)


def gate_result(guard_passed, prefix_passed):
    # Hash equality excluding HEAD is NOT a substitute for the actual guard.
    return {"preflight_passed": bool(guard_passed and prefix_passed),
            "cost_audit": "not_run", "continuation_authorized": False,
            "interval128_executed": False}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def inspect_prefix(d, frozen):
    from sloshing.multiphase.phase_rate_history import Journal, load_archive, read_json
    from sloshing.multiphase.nonlinear_accuracy import NonlinearControls

    folder = d.ROOT / "full_coupling_L0"
    require(not (folder / "FAILED.json").exists(), "Scientific failure marker exists")
    require(not (folder / "COMPLETE.json").exists(), "Prefix is already complete")
    schedule = d.schedules()[0]
    require(schedule.nsteps == 1068 and schedule.sha256 == EXPECTED_SCHEDULE,
            "Unexpected L0 schedule")
    require(schedule.descriptor == frozen["L0_schedule"], "Frozen descriptor differs")
    latest = read_json(folder / "LATEST.json")
    require(latest["checkpoint"] == "checkpoints/000127", "Latest is not checkpoint127")
    arrays, meta = load_archive(folder / latest["checkpoint"])
    expected_identity = {
        "implementation": frozen["implementation"],
        "source_archive_sha256": frozen["source_archive_sha256"],
        "stack": frozen["stack"], "policy_sha256": frozen["policy_sha256"],
        "input_hashes": frozen["input_hashes"],
        "physical_config_fingerprint": frozen["physical_config_fingerprint"],
        "mesh_fingerprint": frozen["historical_input_provenance"]["mesh_fingerprint"],
        "schedule_sha256": EXPECTED_SCHEDULE, "level": "full_CHNS_L0",
        "rank_count": frozen["stack"]["rank_count"],
    }
    require(meta["identity"] == expected_identity, "Checkpoint scientific identity differs")
    require(meta["accepted"] and meta["accepted_steps"] == 127, "Not accepted checkpoint127")
    require(meta["semantic_sha256"] == latest["checkpoint_sha256"], "Latest hash mismatch")
    require(meta["time"] == 6.148107277665776e-10 == schedule.intervals[126].end,
            "Checkpoint time mismatch")
    rows = Journal.read(folder / "scalar_history.jsonl")
    timing = Journal.read(folder / "timing_history.jsonl")
    require([r["step"] for r in rows] == list(range(128)), "Scalar history count/sequence")
    require([r["step"] for r in timing] == list(range(1, 128)), "Timing count/sequence")
    require(meta["journal_prefix"] == {
        "bytes": (folder / "scalar_history.jsonl").stat().st_size,
        "sha256": digest(folder / "scalar_history.jsonl"), "records": 128},
        "Checkpoint does not cover the exact scalar journal")
    for row in rows:
        require(all(row["physical_gates"].values()) and all(row["full_physical_checks"].values()),
                "Historical accepted state fails recorded gates")
    for row, interval in zip(rows[1:], schedule.intervals):
        require((row["time"], row["dt"], row["block"]) ==
                (interval.end, interval.dt, interval.block), "Historical clock mismatch")
        require(row["solver"]["reason"] > 0 and row["solver"]["dt_reductions"] == 0 and
                row["solver"]["actual_controls"] == NonlinearControls().record(),
                "Historical controls/reason mismatch")
    sessions = sorted(folder.glob("sessions/*/session_status.json"))
    require(bool(sessions) and all(read_json(p)["status"] == "complete" for p in sessions),
            "Primary session not closed")
    archives = sorted(folder.glob("fields/*/metadata.json")) + sorted(folder.glob("checkpoints/*/metadata.json"))
    for p in archives:
        load_archive(p.parent, expected_identity)
    require(all((folder / "fields" / f"{i:06d}" / "metadata.json").exists()
                for i in schedule.field_indices if i <= 127), "Scheduled prefix field missing")
    paths = [folder / "scalar_history.jsonl", folder / "timing_history.jsonl", folder / "LATEST.json",
             d.ROOT / "policy/production_frozen.json", d.ROOT / "policy/production_policy.json",
             d.ROOT / "policy/production_source.tar.gz", d.ROOT / "policy/cost_before_full_L0.json",
             d.ROOT / "pilot/status.json", d.ROOT / "STOP.json", Path("STEP3A9_REPORT.md")]
    paths += sessions + archives + [p.parent / "arrays.npz" for p in archives]
    return {"passed": True, "accepted_steps": 127, "expected_steps": 1068,
            "time": meta["time"], "latest_checkpoint": latest["checkpoint"],
            "checkpoint_semantic_sha256": meta["semantic_sha256"],
            "identity": expected_identity, "schedule_descriptor": schedule.descriptor,
            "scalar_history_hash_chain": "passed", "timing_history_hash_chain": "passed",
            "retained_rate_sha256": meta["array_sha256"]["phase_rate"],
            "FAILED_exists": False, "COMPLETE_exists": False,
            "files": bindings(paths), "read_only_diagnostic_after_guard_failure": True}


def run():
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    from sloshing.multiphase.phase_rate_history import atomic_json, read_json
    from sloshing.multiphase.linearized_ch import json_hash

    require(not (ROOT / "preflight.json").exists(), "Refusing to overwrite STEP3A10 preflight")
    # The wrapper mounts sloshing_visualization only, not the repository's .git.
    # It supplies host `git rev-parse HEAD`; read it, never override it.
    actual_head = os.environ.get("STEP3_GIT_COMMIT", "unavailable")
    require(actual_head == BASE, "Unexpected STEP3A10 execution base")
    started = datetime.now(timezone.utc).isoformat()
    guard_error = None
    try:
        d.source_guard()  # Mandatory actual guard. No environment override, no monkeypatch.
    except Exception as exc:
        guard_error = {"type": type(exc).__name__, "message": str(exc)}
    frozen = read_json(d.ROOT / "policy/production_frozen.json")
    actual = d.implementation()
    checks = {
        "implementation": frozen["implementation"] == actual,
        "stack": frozen["stack"] == d.stack(),
        "policy": frozen["policy"] == d.POLICY,
        "policy_sha256": frozen["policy_sha256"] == json_hash(d.POLICY),
        "design_sha256": frozen["design_sha256"] == digest("STEP3A9_DESIGN.md"),
        "schedule": frozen["L0_schedule"] == d.schedules()[0].descriptor,
        "inherited_policy": frozen["inherited_policy"] == read_json("validation_results/step3a4/policy/policy.json"),
        "source_archive_sha256": frozen["source_archive_sha256"] == digest(d.ROOT / "policy/production_source.tar.gz"),
        "audit_proof_sha256": frozen["audit_proof_sha256"] == json_hash(d.audit_proof()),
    }
    with tarfile.open(d.ROOT / "policy/production_source.tar.gz") as tar:
        for member in tar.getmembers():
            if member.isfile():
                path = d.SOURCE / member.name[len("multiphase/"):] if member.name.startswith("multiphase/") else Path(member.name)
                require(path.read_bytes() == tar.extractfile(member).read(), "Archived source bytes differ: " + str(path))
    prefix = inspect_prefix(d, frozen)
    prefix.update(execution_base_HEAD=actual_head, captured_at_UTC=started,
                  source_guard_passed=guard_error is None)
    atomic_json(ROOT / "inherited_prefix_manifest.json", prefix)
    verify_bindings(prefix["files"])
    result = {"execution_base_HEAD": actual_head, "timestamp_UTC": started,
              "source_guard": {"passed": guard_error is None, "error": guard_error},
              "checks": checks, "source_archive_members_byte_identical": True,
              "implementation_differences": {k: {"frozen": frozen["implementation"].get(k), "current": actual.get(k)}
                    for k in sorted(set(frozen["implementation"]) | set(actual))
                    if frozen["implementation"].get(k) != actual.get(k)},
              "current_implementation": actual, "stack": d.stack(),
              "preflight_script_sha256": digest(__file__),
              "prefix_manifest_sha256": digest(ROOT / "inherited_prefix_manifest.json"),
              **gate_result(guard_error is None, prefix["passed"])}
    atomic_json(ROOT / "preflight.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["preflight_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(run())
