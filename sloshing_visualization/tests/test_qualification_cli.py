import importlib.util
from pathlib import Path
import sys

import pytest


def test_standalone_particles_cannot_bypass_full_field_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "path", list(sys.path))
    path = Path(__file__).resolve().parents[1]/"scripts/qualify_animation_dataset.py"
    spec = importlib.util.spec_from_file_location("qualification_cli", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "compare_datasets", lambda *args: {})
    monkeypatch.setattr(module, "fine_short_gate", lambda _: {
        "numerical_gate_passed": False, "problems": ["deliberately rejected field comparison"]})
    monkeypatch.setattr(sys, "argv", [str(path), "--phase", "particles", "--output", str(tmp_path)])
    with pytest.raises(RuntimeError, match="before particles"):
        module.main()
    assert (tmp_path/"full_numerical_gate.json").exists()
    assert not (tmp_path/"runs").exists()
