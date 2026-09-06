import numpy as np
import pytest
from sloshing.visualization.scope import magnify_positions,magnify_surface
from sloshing.visualization.tracers import choose_magnification


def test_display_transforms_are_pure_and_isotropic():
    seed=np.array([[.2,-.3],[-.5,-1.]])
    actual=seed+np.array([[.001,.002],[-.003,.004]])
    original=actual.copy();initial=seed.copy()
    result=magnify_positions(seed,actual,125)
    np.testing.assert_allclose(result-seed,125*(actual-seed))
    np.testing.assert_array_equal(actual,original);np.testing.assert_array_equal(seed,initial)
    eta=np.array([-.003,.001]);before=eta.copy()
    np.testing.assert_array_equal(magnify_surface(eta,10),eta*10)
    np.testing.assert_array_equal(eta,before)
    with pytest.raises(ValueError):magnify_positions(seed,actual,0)


def test_one_factor_from_whole_dataset_not_each_frame():
    d=np.zeros((2,3,2));d[1,:,0]=.001
    m,q=choose_magnification(d,.075)
    assert m==75 and q==.001
