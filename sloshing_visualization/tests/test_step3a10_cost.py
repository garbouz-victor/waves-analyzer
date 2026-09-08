"""Only first-hard-gate/provenance regressions: cost model is NOT RUN."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "step3a10_preflight", Path(__file__).resolve().parents[1] / "scripts/step3a10_preflight.py")
preflight = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preflight)


@pytest.mark.parametrize("guard,prefix", [(False, False), (False, True), (True, False)])
def test_failed_first_gate_never_authorizes_cost_or_pde(guard, prefix):
    result = preflight.gate_result(guard, prefix)
    assert not result["preflight_passed"]
    assert result["cost_audit"] == "not_run"
    assert result["continuation_authorized"] is False
    assert result["interval128_executed"] is False


def test_even_passed_preflight_does_not_substitute_for_cost_authorization():
    result = preflight.gate_result(True, True)
    assert result["preflight_passed"]
    assert not result["continuation_authorized"]
    assert result["cost_audit"] == "not_run"


def test_prefix_binding_detects_same_size_mutation(tmp_path):
    path = tmp_path / "history.jsonl"
    path.write_bytes(b"initial\n")
    records = preflight.bindings([path])
    preflight.verify_bindings(records)
    path.write_bytes(b"changed\n")
    with pytest.raises(ValueError, match="Inherited file changed"):
        preflight.verify_bindings(records)


def test_prefix_binding_detects_size_mismatch(tmp_path):
    path = tmp_path / "checkpoint"
    path.write_bytes(b"unchanged")
    records = preflight.bindings([path])
    records[str(path)]["bytes"] += 1
    with pytest.raises(ValueError, match="Inherited file changed"):
        preflight.verify_bindings(records)


def test_new_files_are_outside_frozen_multiphase_tree():
    assert "/scripts/" in str(Path(spec.origin))
    source = Path(spec.origin).read_text()
    assert "d.source_guard()" in source
    assert "d.execute(" not in source
    assert "d.full_L0(" not in source
    assert "os.environ[" not in source
