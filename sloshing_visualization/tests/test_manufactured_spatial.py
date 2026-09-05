import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver
from sloshing.validation.spatial import convergence_rows, manufactured_run


def test_manufactured_spatial_smoke_in_default_suite():
    # Independent analytic errors on 8x8, with limits set after measuring MMS.
    row = manufactured_run(8)
    assert row["velocity_l2"] < .003
    assert row["velocity_h1"] < .09
    assert row["pressure_l2"] < .006
    assert row["eta_l2"] < .004
    assert row["traction_tangent_l2"] < .02
    assert row["traction_normal_l2"] < .04


def test_production_cannot_enable_manufactured_forcing_via_config():
    solver = SloshingSolver(SimulationConfig(nx=4, nz=4, t_end=0))
    assert solver.integrator.problem.is_production
    for t in (0, .25, 1):
        np.testing.assert_array_equal(solver.integrator.loads(t), 0)
    assert np.count_nonzero(solver.integrator.initial_state().velocity) == 0


@pytest.fixture(scope="module")
def measured_spatial():
    return convergence_rows()


@pytest.mark.validation
def test_observed_manufactured_spatial_orders(measured_spatial):
    # Set after measuring 3.03 / 2.00 / 2.24 / 3.00 on the final pair.
    limits = {"velocity_l2": 2.4, "velocity_h1": 1.7, "pressure_l2": 1.5, "eta_l2": 2.3}
    for row in measured_spatial[1:]:
        for key, limit in limits.items():
            assert row[key+"_order"] > limit


@pytest.mark.validation
def test_mms_spatial_errors_are_not_time_error(measured_spatial):
    row = measured_spatial[-1]
    finer_time = manufactured_run(32, row["dt"]/2)
    for key in ("velocity_l2", "velocity_h1", "pressure_l2", "eta_l2"):
        assert abs(finer_time[key]-row[key])/finer_time[key] < .01
