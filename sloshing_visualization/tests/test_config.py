from dataclasses import replace

import pytest

from sloshing import SimulationConfig


@pytest.mark.parametrize("changes", [{"nu": 0}, {"nu": -1}, {"dt": 0}, {"t_end": -.1},
                                   {"alpha_deg": float("nan")}, {"nx": 2}, {"nx": 8.5},
                                   {"snapshot_dt": .003}, {"t_end": .001}, {"mesh": "unknown"}])
def test_reject_invalid_config(changes):
    with pytest.raises(ValueError):
        replace(SimulationConfig(), **changes)


def test_config_roundtrip_and_stable_run_id(tmp_path):
    config = SimulationConfig()
    path = tmp_path / "config.json"
    path.write_text(config.to_json())
    assert SimulationConfig.from_json(path) == config
    assert replace(config, nu=.1).run_name != config.run_name
