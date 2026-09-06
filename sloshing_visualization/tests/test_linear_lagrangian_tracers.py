import numpy as np
import pytest
from sloshing.visualization.tracers import linear_lagrangian_displacement


def test_zero_and_uniform_time_linear_velocity():
    t=np.linspace(0,1,11)
    v=np.zeros((11,3,2))
    assert not linear_lagrangian_displacement(t,v).any()
    v[:]=t[:,None,None]*np.array([2.,-3.])
    exact=t[:,None,None]**2*np.ones((1,3,1))*np.array([1.,-1.5])
    np.testing.assert_allclose(linear_lagrangian_displacement(t,v),exact,atol=3e-16)
    # Spatially uniform advection has this identical integral at displaced X.


def test_smooth_integral_second_order_at_fixed_reference_points():
    errors=[]
    x0=np.array([-.5,.25])
    for n in (20,40,80):
        t=np.linspace(0,1,n+1)
        v=np.stack((np.exp(t[:,None])*x0[None],np.zeros((n+1,2))),axis=-1)
        xi=linear_lagrangian_displacement(t,v)
        errors.append(np.max(abs(xi[-1,:,0]-(np.e-1)*x0)))
    assert min(np.log2(np.array(errors[:-1])/errors[1:]))>1.99


def test_spatially_varying_paths_differ_only_at_higher_order():
    # v=epsilon*x. First-order displacement epsilon*x0*t;
    # fully advected position x0*exp(epsilon*t). Difference ~epsilon².
    errors=[]
    for eps in (.04,.02,.01):
        t=np.linspace(0,1,11)
        v=np.zeros((11,1,2));v[:,:,0]=eps*.5
        linear=.5+linear_lagrangian_displacement(t,v)[-1,0,0]
        errors.append(abs(.5*np.exp(eps)-linear))
    assert min(np.log2(np.array(errors[:-1])/errors[1:]))>1.99


def test_invalid_time_grid_rejected():
    with pytest.raises(ValueError):linear_lagrangian_displacement([0,0],np.zeros((2,1,2)))


def test_uniform_linear_tracer_matches_actual_rk_pathline():
    from sloshing.validation.particles import advect
    class UniformHistory:
        times=np.linspace(0,1,101)
        def valid(self,t,p):return np.ones(p.shape[1],dtype=bool)
        def velocity(self,t,p):return np.tile(np.array([.01*t,-.02*t])[:,None],(1,p.shape[1]))
    hist=UniformHistory();seeds=np.array([[.2,-.3],[-.5,-1.]])
    v=np.array([hist.velocity(t,seeds.T).T for t in hist.times])
    linear=seeds+linear_lagrangian_displacement(hist.times,v)
    path=advect(hist,seeds,max_step=.005)["paths"]
    np.testing.assert_allclose(path,linear,atol=2e-15,rtol=0)
