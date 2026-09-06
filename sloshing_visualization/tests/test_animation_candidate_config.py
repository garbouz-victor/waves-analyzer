from dataclasses import asdict
from pathlib import Path

from sloshing.config import SimulationConfig


CONFIGS = Path(__file__).resolve().parents[1]/"configs"


def test_animation_candidates_and_historical_config():
    c = SimulationConfig.from_json(CONFIGS/"animation_candidate.json")
    f = SimulationConfig.from_json(CONFIGS/"animation_candidate_fine.json")
    assert c.alpha_deg == .02 and c.nu == .01 and c.t_end == 5
    assert c.integrator == "sdirk2" and c.dt == .00125 and c.snapshot_dt == .005
    assert c.resolution == (80, 140) and f.resolution == (120, 210)
    assert {k:v for k,v in asdict(c).items() if k != "mesh"} == {k:v for k,v in asdict(f).items() if k != "mesh"}
    old = SimulationConfig.from_json(CONFIGS/"validation_alpha_0_2.json")
    assert old.alpha_deg == .2 and old.integrator == "midpoint" and old.t_end == .5
    assert SimulationConfig().alpha_deg == 2 and SimulationConfig().integrator == "midpoint"
    assert not (CONFIGS/"linear_safe.json").exists()
