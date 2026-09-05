import numpy as np


def test_volume_is_preserved_including_endpoint_basis_functions(short_run):
    solver, snapshots = short_run
    f = solver.fem
    assert abs(f.surface_weights.sum() - 2 * solver.config.a) < 1e-13
    for state, row in snapshots:
        assert abs(row["surface_integral"]) < 1e-12
        # Independent exact Simpson integration on every quadratic surface edge.
        eta, x = state.eta, f.surface_x
        integral = np.sum((x[2::2] - x[:-2:2]) / 6 * (eta[:-2:2] + 4*eta[1::2] + eta[2::2]))
        assert abs(integral) < 1e-12
