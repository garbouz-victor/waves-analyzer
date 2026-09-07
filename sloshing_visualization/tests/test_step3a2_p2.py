import numpy as np
import pytest
from sloshing.multiphase.resolution import REFERENCE_POINTS, polynomial_bounds, active_transition_cells
from sloshing.multiphase.p2_certification import polynomial_ranges, certified_normal_widths, coefficients, evaluate,physical_gradients


def measure(fun):
    x,y=REFERENCE_POINTS[:6].T
    return polynomial_ranges(np.asarray(fun(x,y)+np.zeros(6))[None,:])


@pytest.mark.parametrize("fun,lo,hi,mintype,maxtype",[
    (lambda x,y:2.,2.,2.,"vertex","vertex"),
    (lambda x,y:2*x-3*y+1,-2.,3.,"vertex","vertex"),
    (lambda x,y:(x-.2)**2+(y-.3)**2,0.,.73,"interior","vertex"),
    (lambda x,y:-(x-.2)**2-(y-.3)**2,-.73,0.,"vertex","interior"),
    (lambda x,y:1-(x-.5)**2-y,-.25,1.,"vertex","edge"),
    (lambda x,y:(x-.5)**2+y,0.,1.25,"edge","vertex"),
    (lambda x,y:(x-.2)**2-(y-.2)**2,-.6,.6,"vertex","vertex"),
    (lambda x,y:(x+y-.5)**2,0.,.25,"edge","vertex"),
])
def test_exact_extrema_and_witnesses(fun,lo,hi,mintype,maxtype):
    r=measure(fun)
    assert r["exact"][0]
    assert r["min"][0]==pytest.approx(lo,abs=2e-12)
    assert r["max"][0]==pytest.approx(hi,abs=2e-12)
    assert r["min_type"][0]==mintype and r["max_type"][0]==maxtype
    for side in ("min","max"):
        x,y=r[side+"_point"][0]
        assert x>=0 and y>=0 and x+y<=1+1e-15
        assert fun(x,y)==pytest.approx(r[side+"_value"][0],abs=1e-12)


def test_nearly_singular_hessian_uses_labelled_outer_bound():
    r=measure(lambda x,y:(x-.25)**2+1e-14*(y-.25)**2)
    assert not r["exact"][0]
    assert r["status"][0]=="certified_conservative_fallback"
    assert r["min"][0]<=0


@pytest.mark.parametrize("fun",[
    lambda x,y:10*(x+y)-1,
    lambda x,y:1.2-1.6*x*(1-x-y),
])
def test_true_transition_near_vertex_or_edge_with_pure_centroid(fun):
    assert fun(1/3,1/3)>.9
    r=measure(fun)
    assert r["min"][0]<=.9 and r["max"][0]>=-.9


def test_bernstein_false_positive_removed():
    x,y=REFERENCE_POINTS[:6].T
    values=(1.+(x-.5)**2)[None,:]
    assert polynomial_bounds(values,2)[0][0]<.9
    assert not active_transition_cells(values,2)[0]
    assert measure(lambda x,y:1.+(x-.5)**2)["min"][0]>.9


def test_random_polynomials_never_escape_exact_bounds():
    rng=np.random.default_rng(3291)
    values=rng.normal(size=(350,6))
    co=coefficients(values)
    assert evaluate(co[:,None,:],REFERENCE_POINTS[None,:6,:])==pytest.approx(values,abs=1e-12)
    r=polynomial_ranges(values)
    grid=np.array([(i/100,j/100) for i in range(101) for j in range(101-i)])
    sampled=evaluate(co[:,None,:],grid[None,:,:])
    assert np.all(sampled.min(axis=1)>=r["min"]-1e-12)
    assert np.all(sampled.max(axis=1)<=r["max"]+1e-12)


def test_p1_is_exact_at_vertices():
    r=polynomial_ranges(np.array([[1.,-1.,4.]]),1)
    assert r["min"][0]==-1. and r["max"][0]==4.


TRI=np.array([[[0.,0.],[1.,0.],[0.,1.]]])


def test_constant_gradient_exact_projected_width():
    g=np.tile([[[.4,.9]]],(1,3,1))
    r=certified_normal_widths(TRI,g)
    expected=np.ptp(TRI[0]@([.4,.9]/np.linalg.norm([.4,.9])))
    assert not r["diameter_fallback"][0]
    assert r["width"][0]==pytest.approx(expected,abs=3e-12)


@pytest.mark.parametrize("angle",[.2,1.5,2.9,5.9])
def test_width_rotation_covariance(angle):
    R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    g=np.array([[[.1,1.],[.2,1.],[.3,1.]]])
    assert certified_normal_widths(TRI@R.T,g@R.T)["width"]==pytest.approx(certified_normal_widths(TRI,g)["width"],abs=2e-12)


@pytest.mark.parametrize("g",[
    [[[1.,0.],[0.,1.],[-1.,-1.]]],
    [[[1.,0.],[-1.,0.],[2.,0.]]],
    [[[0.,0.],[0.,0.],[0.,0.]]],
    [[[1.,1e-13],[-1.,1e-13],[2.,1e-13]]],
])
def test_origin_or_near_origin_uses_diameter(g):
    r=certified_normal_widths(TRI,np.array(g))
    assert r["diameter_fallback"][0]
    assert r["width"][0]==pytest.approx(np.sqrt(2))


def test_random_gradient_hulls_never_underestimate_width():
    rng=np.random.default_rng(8012)
    tri=rng.normal(size=(250,3,2))
    g=rng.normal(size=(250,3,2))
    # Also exercise many narrow cones far from the origin.
    g[125:]=np.array([1.,2.])+g[125:]*.02
    r=certified_normal_widths(tri,g)
    weights=rng.dirichlet([1,1,1],4000)
    samples=np.einsum("pv,cvi->cpi",weights,g)
    samples/=np.linalg.norm(samples,axis=-1,keepdims=True)
    widths=np.ptp(np.einsum("cvi,cpi->cpv",tri,samples),axis=-1)
    assert np.all(widths<=r["width"][:,None]+1e-12)
    assert (~r["diameter_fallback"])[125:].all()


def test_physical_gradient_from_coefficients_is_exact_for_constant_and_linear():
    tri=np.array([[[.1,.2],[1.3,.3],[.2,1.4]]])
    ref=REFERENCE_POINTS[:6]
    physical=tri[:,0,None,:]+np.einsum("pi,cij->cpj",ref,tri[:,1:]-tri[:,0,None,:])
    for degree in (1,2):
        n=3 if degree==1 else 6
        constant=physical_gradients(tri,np.full((1,n),.345),REFERENCE_POINTS,degree)
        assert np.max(abs(constant))==0.
        values=3*physical[:,:n,0]-2*physical[:,:n,1]+.4
        gradient=physical_gradients(tri,values,REFERENCE_POINTS,degree)
        assert gradient==pytest.approx(np.broadcast_to([3.,-2.],gradient.shape),abs=1e-13)
