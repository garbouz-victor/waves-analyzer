"""Equal zero-contour angle does not guarantee the full diffuse natural BC."""
import sympy as sp


def test_curved_tanh_wall_variation_derived_symbolically():
    x,z,R,eps,sigma=sp.symbols("x z R eps sigma",positive=True)
    radius=sp.sqrt(x*x+(z+R/2)**2)  # theta_geometry=60 degrees
    phi=sp.tanh((R-radius)/(sp.sqrt(2)*eps))
    lam=3*sigma/(2*sp.sqrt(2))
    wall=-3*sigma*(1-phi**2)/8  # f_wall' at theta_target=60 degrees
    actual=(lam*eps*(-sp.diff(phi,z))+wall).subs(z,0)
    r0=radius.subs(z,0);phi0=phi.subs(z,0)
    expected=3*sigma/8*(R/r0-1)*(1-phi0**2)
    assert sp.simplify(actual-expected)==0
    assert sp.simplify(expected.subs(x,sp.sqrt(3)*R/2))==0  # at phi=0 only
    assert abs(float(expected.subs({R:.3,eps:.025,sigma:.1,x:.28})))>1e-4
