"""Authorized external numerical-identity guard; historical guard is untouched."""
import hashlib
import json
import os
from pathlib import Path
import tarfile
from datetime import datetime, timezone

ROOT = Path("validation_results/step3a10")
OLD = Path("validation_results/step3a9")
TRAJECTORY = OLD / "full_coupling_L0"
ADMIN_VERSION = "step3a10-numerical-identity-with-separate-session-lineage-v1"
MANIFEST = ROOT / "inherited_prefix_manifest.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def semantic(value):
    # Match the existing canonical JSON helper, not a new numerical identity.
    from sloshing.multiphase.linearized_ch import json_hash
    return json_hash(value)


def numerical_implementation(implementation):
    if "execution_base_HEAD" not in implementation:
        raise ValueError("Missing explicit execution lineage")
    return {k: v for k, v in implementation.items() if k != "execution_base_HEAD"}


def require_numerical_match(frozen, current):
    if numerical_implementation(frozen) != numerical_implementation(current):
        raise ValueError("Numerical implementation identity changed")


def numerical_guard():
    """Same frozen checks; ONLY execution_base_HEAD is recorded separately.

    No substitution into the old module and no environment assignment. This
    callback can be passed through the frozen runner's public source_guard API.
    """
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    f = read(OLD / "policy/production_frozen.json")
    current = d.implementation()
    require_numerical_match(f["implementation"], current)
    if current["execution_base_HEAD"] == "unavailable":
        raise ValueError("Actual session HEAD must be supplied by unchanged wrapper")
    checks = {
        "stack": f["stack"] == d.stack(),
        "policy": f["policy"] == d.POLICY,
        "policy_sha256": f["policy_sha256"] == semantic(d.POLICY),
        "design": f["design_sha256"] == sha("STEP3A9_DESIGN.md"),
        "schedule": f["L0_schedule"] == d.schedules()[0].descriptor,
        "inherited_policy": f["inherited_policy"] == read("validation_results/step3a4/policy/policy.json"),
        "source_archive": f["source_archive_sha256"] == sha(OLD / "policy/production_source.tar.gz"),
        "audit_proof": f["audit_proof_sha256"] == semantic(d.audit_proof()),
    }
    if not all(checks.values()):
        raise ValueError("Numerical guard: " + str(checks))
    return f


def lineage(frozen):
    actual = os.environ.get("STEP3_GIT_COMMIT", "unavailable")
    if actual == "unavailable":
        raise ValueError("Missing actual wrapper HEAD")
    return {"historical_execution_base_HEAD": frozen["implementation"]["execution_base_HEAD"],
            "current_session_HEAD": actual, "administrative_guard_version": ADMIN_VERSION,
            "numerical_identity_sha256": semantic(numerical_implementation(frozen["implementation"])),
            "historical_checkpoint_identity_role": "immutable trajectory origin, NOT current session HEAD"}


def prefix_guard(*, after_continuation=False):
    """Before: exact files. After: byte-preserved journal prefixes/old archives."""
    manifest = read(MANIFEST)
    for name, record in manifest["files"].items():
        p = Path(name)
        if after_continuation and p == TRAJECTORY / "LATEST.json":
            continue  # Only the frozen runner may update this accepted pointer.
        if after_continuation and p in (TRAJECTORY / "scalar_history.jsonl", TRAJECTORY / "timing_history.jsonl"):
            with p.open("rb") as stream:
                data = stream.read(record["bytes"])
            if hashlib.sha256(data).hexdigest() != record["sha256"]:
                raise ValueError("Historical journal bytes changed: " + name)
        elif sha(p) != record["sha256"] or p.stat().st_size != record["bytes"]:
            raise ValueError("Historical immutable file changed: " + name)
    if not after_continuation:
        from sloshing.multiphase.benchmarks import divergence_transfer as d
        from step3a10_preflight import inspect_prefix
        result = inspect_prefix(d, read(OLD / "policy/production_frozen.json"))
        if not result["passed"]:
            raise ValueError("Prefix127 integrity failed")
    return manifest


def require_authorization(root=ROOT):
    """Pure file binding gate. Reject BEFORE importing the PDE execution path."""
    root = Path(root)
    auth = read(root / "continuation_authorization.json")
    if auth.get("authorized") is not True:
        raise ValueError("Continuation is not authorized")
    for relative, key in (("cost_audit/policy.json", "cost_policy_sha256"),
                          ("cost_audit/decision.json", "cost_decision_sha256"),
                          ("inherited_prefix_manifest.json", "prefix_manifest_sha256")):
        if sha(root / relative) != auth[key]:
            raise ValueError("Authorization binding mismatch: " + key)
    policy, decision = read(root / "cost_audit/policy.json"), read(root / "cost_audit/decision.json")
    if decision.get("authorized") is not True or decision["policy_sha256"] != auth["cost_policy_sha256"]:
        raise ValueError("Decision/policy is not authorized")
    if (decision["prefix_manifest_sha256"] != auth["prefix_manifest_sha256"] or
            decision["schedule_sha256"] != auth["schedule_sha256"]):
        raise ValueError("Decision/prefix/schedule binding mismatch")
    if (decision["forecast"]["projected_complete_L0_s"] > policy["full_L0_cap_s"] or
            decision["forecast"]["projected_total_scientific_s"] > policy["total_cap_s"] or
            not all(row["passed"] for row in decision["backtests"])):
        raise ValueError("Failed frozen cost decision")
    for path, digest in auth["administrative_files_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Authorized administrative source changed: " + path)
    return auth


def record_guard():
    from sloshing.multiphase.phase_rate_history import atomic_json
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    f = numerical_guard()  # No cost analysis precedes this call.
    prefix_guard()
    members = 0
    with tarfile.open(OLD / "policy/production_source.tar.gz") as tar:
        for member in tar.getmembers():
            if member.isfile():
                p = d.SOURCE / member.name[len("multiphase/"):] if member.name.startswith("multiphase/") else Path(member.name)
                if p.read_bytes() != tar.extractfile(member).read():
                    raise ValueError("Frozen archive member changed: " + str(p))
                members += 1
    result = {"passed": True, "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
              **lineage(f), "source_archive_members_byte_identical": members,
              "prefix_manifest_sha256": sha(MANIFEST), "prefix_integrity": "passed",
              "source_archive_sha256": f["source_archive_sha256"], "policy_sha256": f["policy_sha256"],
              "schedule_sha256": f["L0_schedule"]["schedule_sha256"],
              "guard_source_sha256": sha(__file__), "old_guard_modified": False,
              "environment_modified": False, "checkpoint127_modified": False}
    out = ROOT / "administrative_guard/initial_PASS.json"
    if out.exists():
        raise ValueError("Refusing to overwrite administrative guard evidence")
    atomic_json(out, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    record_guard()
