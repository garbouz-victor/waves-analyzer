"""Pure I/O/policy tests only: no FEM or trajectory solver construction."""
import ast
import copy
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from step3a11_io_attribution import (InstrumentedJournal, Journal, PHASES, quantiles,
    wilson, operation_classification, context_classification, check_bindings, sha,
    instrumented_atomic, atomic_json, safe_cleanup)


def payloads():
    return [{"step": i, "time": i*4.657657028534679e-12, "solver": {"reason": 2, "trace": [1e-4, 1e-11]},
             "nested": {"values": [True, None, "μ", -0., 1.234567890123456e-15]}, "energy": -1e-7+i*1e-13}
            for i in range(10)]


def test_ten_byte_identical_append_chains(tmp_path):
    a, b = tmp_path/"a.jsonl", tmp_path/"b.jsonl"
    a.write_bytes(b""); b.write_bytes(b"")
    original = Journal(a); measured = InstrumentedJournal(b)
    for i, row in enumerate(payloads()):
        before = copy.deepcopy(row)
        original.append(row); _, rec = measured.append(row, probes=False, repetition=i)
        assert row == before
        assert a.read_bytes() == b.read_bytes()
        assert Journal.read(a) == Journal.read(b)
        assert rec["file_size_after"]-rec["file_size_before"] == rec["line_bytes"]
        assert sum(rec[k+"_s"] for k in PHASES) <= rec["total_append_s"]+1e-9
        assert all(rec[k+"_s"] >= 0 for k in PHASES)


def sample(phase="fsync", slow=8., other=.001):
    row = {k+"_s": other for k in PHASES}
    row[phase+"_s"] = slow
    row["total_append_s"] = sum(row.values())+.001
    row["non_fsync_s"] = row["total_append_s"]-row["fsync_s"]
    return row


def test_fsync_dominated():
    result = operation_classification([sample() for _ in range(5)])
    assert result["classification"] == "FSYNC_DOMINATED"
    assert all(result["checks"].values())


@pytest.mark.parametrize("phase", ["write", "flush", "open", "close", "json_encode", "record_hash", "copy_prepare"])
def test_other_phase_is_not_fsync(phase):
    assert operation_classification([sample(phase) for _ in range(10)])["classification"] == "FSYNC_NOT_YET_ISOLATED"


@pytest.mark.parametrize("n", [0, 1, 4])
def test_fewer_than_five_stalls(n):
    assert operation_classification([sample() for _ in range(n)])["classification"] == "FSYNC_NOT_YET_ISOLATED"


def test_non_fsync_threshold_is_not_relaxed():
    assert not operation_classification([sample(other=.02) for _ in range(8)])["checks"]["non_fsync_p95"]


def metrics(values, fractions):
    return {name:{"p95":x, "stall_fraction":f} for name,x,f in zip(
        ("bind", "host", "container_tmp", "docker_volume"), values, fractions)}


def test_bind_context():
    assert context_classification(metrics([8., .01, .01, .01], [.2, 0., 0., 0.]), True)["classification"] == "BIND_MOUNT_SPECIFIC_EVIDENCE"


def test_host_context():
    assert context_classification(metrics([8., 7., .01, .01], [.2, .2, 0., 0.]), True)["classification"] == "HOST_FILESYSTEM_DURABILITY_EVIDENCE"


def test_docker_context():
    assert context_classification(metrics([8., .01, 7., 8.], [.2, 0., .2, .2]), True)["classification"] == "DOCKER_STORAGE_PATH_EVIDENCE"


@pytest.mark.parametrize("operation_pass", [False, True])
def test_ambiguous_context(operation_pass):
    assert context_classification(metrics([8., 7., 8., 9.], [.2]*4), operation_pass)["classification"] == "STORAGE_PATH_UNRESOLVED"


def test_operation_precedes_path():
    assert context_classification(metrics([8., .01, .01, .01], [.2, 0., 0., 0.]), False)["classification"] == "STORAGE_PATH_UNRESOLVED"


def test_wilson_and_percentiles():
    low, high = wilson(22, 109)
    assert .12 < low < 22/109 < high < .30
    assert wilson(0, 120)[0] == 0
    assert quantiles([1., 2., 3., 4., 5.])["p95"] == pytest.approx(4.8)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_nonfinite_metrics_rejected(bad):
    with pytest.raises(ValueError): quantiles([1., bad])


def test_atomic_byte_equivalence(tmp_path):
    a,b = tmp_path/"a.json", tmp_path/"b.json"
    atomic_json(a, payloads()[3]); result = instrumented_atomic(b, payloads()[3], probes=False)
    assert a.read_bytes() == b.read_bytes()
    assert result["fsync_s"] >= 0


def test_binding_detects_mutation(tmp_path):
    p = tmp_path/"immutable"; p.write_bytes(b"unchanged")
    bindings = {str(p):{"sha256":sha(p), "bytes":p.stat().st_size}}
    check_bindings(bindings); p.write_bytes(b"changed")
    with pytest.raises(ValueError): check_bindings(bindings)


def test_cleanup_refuses_unmarked_directory(tmp_path):
    p = tmp_path/"step3a11-test"; p.mkdir()
    with pytest.raises(ValueError): safe_cleanup(p, tmp_path)
    assert p.exists()


def test_no_numerical_constructors_or_steps():
    scripts = Path(__file__).resolve().parents[1]/"scripts"
    forbidden = {"CHNSSolver", "CHNSPhaseRateBE", "DivergenceTrajectory", "RateTrajectory", "advance", "run_until", "solve"}
    for name in ("step3a11_io_attribution.py", "step3a11_host_io.py", "step3a11_report.py"):
        tree = ast.parse((scripts/name).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                called = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ""
                assert called not in forbidden


def test_current_prefix_still_127_and_old_authorization_false():
    root = Path(__file__).resolve().parents[1]/"validation_results"
    prefix = root/"step3a9/full_coupling_L0"
    assert [r["step"] for r in Journal.read(prefix/"scalar_history.jsonl")] == list(range(128))
    assert json.loads((prefix/"LATEST.json").read_text())["checkpoint"] == "checkpoints/000127"
    assert len(list(prefix.glob("sessions/*/session_status.json"))) == 2
    for name in ("FAILED.json", "COMPLETE.json", "fields/000128", "checkpoints/000128"):
        assert not (prefix/name).exists()
    assert json.loads((root/"step3a10/continuation_authorization.json").read_text())["authorized"] is False
