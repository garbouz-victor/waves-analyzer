import numpy as np

from sloshing.config import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.validation.qualification.fem_evaluation import FEMGeometry
from sloshing.validation.qualification.signals import horizontal_quadrature, modal_events


def test_modal_event_estimator_on_independent_damped_signal():
    t = np.linspace(0, 5, 1001)
    period, decay = 1.6, .1
    omega = 2*np.pi/period
    result = modal_events(t, np.exp(-decay*t)*np.cos(omega*t))
    assert abs(result["estimated_period_s"]-period) < 1e-5
    assert abs(result["first_zero_crossing_s"]-period/4) < 1e-5
    expected_min = (np.pi-np.arctan(decay/omega))/omega
    assert abs(result["first_opposite_extremum"]["time"]-expected_min) < 3e-5
    for r in result["same_sign_damping"]:
        assert abs(r["effective_decay_rate_per_s"]-decay) < 1e-5


def test_horizontal_fem_line_quadrature_is_exact_for_p2_squared():
    f = FEMSystem(SimulationConfig(nx=4, nz=6))
    geometry = FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs)
    x, z = f.scalar.doflocs
    coeff = x*x+z*x
    q, weights = horizontal_quadrature(geometry, -.13)
    values = geometry.evaluate(coeff, q)
    exact = 2/5+.13**2*2/3  # Integral (x^2-.13*x)^2 on [-1,1].
    np.testing.assert_allclose(weights@values**2, exact, atol=1e-12)
