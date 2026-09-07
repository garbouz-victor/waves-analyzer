import numpy as np
import pytest
from sloshing.multiphase.free_energy import bulk_energy,bulk_derivative,wall_energy,wall_derivative
from sloshing.multiphase.be_energy_work import bulk_remainder,wall_remainder


def test_exact_quartic_BE_remainder_identity():
    rng=np.random.default_rng(2201)
    a,b=rng.uniform(-1.2,1.2,(2,1000))
    actual=bulk_derivative(b,.1,.025)*(b-a)-(bulk_energy(b,.1,.025)-bulk_energy(a,.1,.025))
    assert bulk_remainder(a,b,.1,.025)==pytest.approx(actual,abs=1e-13)


def test_exact_cubic_wall_BE_remainder_identity():
    rng=np.random.default_rng(2202)
    a,b=rng.uniform(-1.2,1.2,(2,1000))
    actual=wall_derivative(b,.1,60)*(b-a)-(wall_energy(b,.1,60)-wall_energy(a,.1,60))
    assert wall_remainder(a,b,.1,60)==pytest.approx(actual,abs=1e-14)


def test_nonconvex_remainder_is_not_a_nonnegative_dissipation():
    assert bulk_remainder(0.,.01,.1,.025)<0
    assert wall_remainder(-1.,-.99,.1,60)<0
