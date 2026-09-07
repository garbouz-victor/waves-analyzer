"""Contour roots, connectivity and contact angle independently of a solver."""
import numpy as np
import pytest
from sloshing.multiphase.interface import quadratic_edge_roots, connected_polylines, fit_circle
from sloshing.multiphase.contact_line import all_wall_crossings, sessile_apparent_angle


def test_quadratic_all_edge_crossings():
    f=lambda s:(s-.2)*(s-.8)
    assert quadratic_edge_roots(f(0),f(.5),f(1)) == pytest.approx([.2,.8])
    assert quadratic_edge_roots(-1,0,1)==pytest.approx([.5])
    assert quadratic_edge_roots(-.0010357373660812043,-.00025903189016354385,
                                3.1401849173675494e-15)==pytest.approx([1.])


def test_all_components_and_multiple_wall_crossings():
    segments=np.array([[[0,0],[1,1]],[[1,1],[0,2]],[[0,3],[1,4]],[[1,4],[0,5]]],float)
    lines=connected_polylines(segments)
    assert len(lines)==2
    result=all_wall_crossings(lines,0.,axis=0)
    assert sorted(r["coordinate"] for r in result)==[0.,2.,3.,5.]


@pytest.mark.parametrize("angle",[60.,90.,120.])
def test_apparent_angle_fit_window_circle(angle):
    radius=.3
    center=-radius*np.cos(np.deg2rad(angle))
    t=np.linspace(0,2*np.pi,800)
    points=np.column_stack([radius*np.cos(t),center+radius*np.sin(t)])
    points=points[points[:,1]>=0]
    for window in ((2,6),(2,10),(3,8)):
        result=sessile_apparent_angle(points,0.,.015,window)
        assert result["theta_deg"]==pytest.approx(angle,abs=1e-10)
        assert result["radial_rms"]<1e-12
