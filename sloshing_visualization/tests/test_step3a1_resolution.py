import numpy as np
import pytest
from sloshing.multiphase.resolution import (REFERENCE_POINTS, polynomial_bounds,
    active_transition_cells,normal_widths,local_resolution)


@pytest.mark.parametrize("values",[[-1.,4.,4.],[-.9,3.,3.],[4.,-1.,4.]])
def test_p1_band_crossing_is_not_lost_with_pure_centroid(values):
    assert np.mean(values)>.9  # old centroid test misses it
    assert active_transition_cells(np.array([values]),1)[0]
    lo,hi=polynomial_bounds(np.array([values]),1)
    assert lo[0]==min(values) and hi[0]==max(values)


def test_p2_edge_transition_with_pure_centroid_and_vertices():
    x,z=REFERENCE_POINTS.T
    phi=1.2-1.6*x*(1-x-z)
    assert phi[-1]>.9 and np.min(phi[:3])>.9
    assert phi[3]<.9
    assert active_transition_cells(phi[None,:6],2)[0]


def test_p2_curved_zero_and_bernstein_convex_hull():
    x,z=REFERENCE_POINTS[:6].T
    phi=50*((x-.1)**2+(z-.1)**2-.08**2)
    assert active_transition_cells(phi[None,:],2)[0]
    lo,hi=polynomial_bounds(phi[None,:],2)
    rng=np.random.default_rng(11)
    points=rng.dirichlet([1.,1.,1.],1000)
    actual=50*((points[:,1]-.1)**2+(points[:,2]-.1)**2-.08**2)
    assert lo[0]<=actual.min() and hi[0]>=actual.max()


def test_normal_width_uses_all_gradients_and_refines_linearly():
    tri=np.array([[[0.,0.],[1.,0.],[0.,1.]]])
    gradients=np.ones((1,7,2));gradients[0,0]=[-1.,1.]
    width=normal_widths(tri,gradients)
    assert width[0]==pytest.approx(np.sqrt(2))
    assert normal_widths(tri/2,gradients)[0]==pytest.approx(width[0]/2)
    values=np.zeros((1,6))
    _,_,counts=local_resolution(tri,values,gradients,2,.1)
    _,_,fine=local_resolution(tri/2,values,gradients,2,.1)
    assert fine[0]==pytest.approx(2*counts[0])


@pytest.mark.parametrize("angle",[.17,.5,1.2,2.5])
def test_rotation_of_physical_mesh_and_interface_preserves_resolution(angle):
    rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    tri=np.array([[[0.,0.],[1.,0.],[0.,1.]]])
    grad=np.tile(np.array([[[.7,-.2]]]),(1,7,1))
    assert normal_widths(tri@rotation.T,grad@rotation.T)==pytest.approx(normal_widths(tri,grad))


def test_constant_phase_in_transition_not_silently_dropped():
    tri=np.array([[[0.,0.],[1.,0.],[0.,1.]]])
    indices,width,_=local_resolution(tri,np.zeros((1,6)),np.zeros((1,7,2)),2,.1)
    assert indices.tolist()==[0] and width[0]==pytest.approx(np.sqrt(2))


def test_interface_rotation_relative_to_fixed_triangle_is_continuous():
    tri=np.array([[[0.,0.],[1.,0.],[0.,1.]]])
    widths=[]
    for angle in np.linspace(0.,2*np.pi,361):
        grad=np.tile([[[np.cos(angle),np.sin(angle)]]],(1,7,1))
        widths.append(normal_widths(tri,grad)[0])
    # Real projected-width anisotropy is allowed; discontinuous centroid on/off
    # artifacts are not. Support width is Lipschitz with the cell diameter.
    assert np.max(np.abs(np.diff(widths)))<=np.sqrt(2)*np.pi/180+1e-12
    assert max(widths)/min(widths)<=2+1e-12
