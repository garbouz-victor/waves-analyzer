import pytest

from sloshing import SimulationConfig
from sloshing.validation.physical import run_history


def test_interrupted_cache_is_neither_reused_nor_overwritten(tmp_path):
    path = tmp_path/"interrupted.h5"
    content = b"truncated HDF5 stand-in"
    path.write_bytes(content)
    with pytest.raises(RuntimeError, match="Unreadable validation cache"):
        run_history(None, SimulationConfig(), path, None, None)
    assert path.read_bytes() == content
