import numpy as np

from sloshing.validation.stiffness import amplification, stiff_fem_comparison


def test_fast_constrained_fem_mode_distinguishes_integrators():
    rows = stiff_fem_comparison()
    for row in rows:
        expected = amplification(row["integrator"], row["lambda_dt"])**row["step"]
        np.testing.assert_allclose(row["amplitude"], expected, atol=1e-10)
        assert row["weak_constraint_max"] < 1e-10
    mid, sdirk = rows[4], rows[9]
    assert mid["amplitude"] < -.8
    assert sdirk["energy_norm_error"] < 1e-6
    assert mid["energy_norm_error"] > .8
