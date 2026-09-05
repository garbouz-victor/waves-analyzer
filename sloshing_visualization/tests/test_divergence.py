import numpy as np

from sloshing import SimulationConfig, SloshingSolver


def test_weak_and_strong_divergence_are_reported_separately(short_run):
    _, snapshots = short_run
    for _, row in snapshots:
        assert row["weak_divergence_l2"] < 1e-11
        assert np.isfinite(row["divergence_l2"])
    assert snapshots[-1][1]["divergence_l2"] > 1e-9  # P2/P1 is not pointwise solenoidal.


def test_actual_divergence_decreases_on_refinement():
    norms = []
    for nx, nz in [(8, 16), (16, 32), (24, 48)]:
        solver = SloshingSolver(SimulationConfig(nx=nx, nz=nz, t_end=.05))
        norms.append(list(solver.snapshots())[-1][1]["divergence_l2"])
    assert norms[1] < norms[0]
    assert norms[2] < norms[1]
