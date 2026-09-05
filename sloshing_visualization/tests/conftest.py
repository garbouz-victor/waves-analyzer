import pytest

from sloshing import SimulationConfig, SloshingSolver


@pytest.fixture(scope="session")
def short_run():
    solver = SloshingSolver(SimulationConfig(nx=16, nz=32, t_end=.1, snapshot_dt=.025,
                                           visualization_nx=9, visualization_nz=13))
    return solver, list(solver.snapshots())
