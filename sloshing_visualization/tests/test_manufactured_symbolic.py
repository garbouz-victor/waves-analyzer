import sympy as sp

from sloshing.validation.manufactured import expressions


def test_symbolic_incompressibility_and_wall_no_slip():
    e = expressions()
    v, x, z = e["velocity"], e["x"], e["z"]
    assert sp.simplify(v[0].diff(x)+v[1].diff(z)) == 0
    for substitution in ({x: -1}, {x: 1}, {z: -1}):
        assert v.subs(substitution).applyfunc(sp.simplify) == sp.zeros(2, 1)


def test_symbolic_kinematics_and_pinned_endpoints():
    e = expressions()
    assert sp.simplify(e["eta"].diff(e["t"])-e["velocity"][1].subs(e["z"], 0)) == 0
    assert e["eta"].subs(e["x"], -1) == e["eta"].subs(e["x"], 1) == 0


def test_symbolic_momentum_including_stress_divergence():
    e = expressions()
    x, z, t, v, q = (e[key] for key in ("x", "z", "t", "velocity", "pressure"))
    residual = v.diff(t)+sp.Matrix([q.diff(x), q.diff(z)])-e["nu"]*(v.diff(x, 2)+v.diff(z, 2))-e["force"]
    assert residual.applyfunc(sp.simplify) == sp.zeros(2, 1)
    div_stress = sp.Matrix([e["stress"][i, 0].diff(x)+e["stress"][i, 1].diff(z) for i in range(2)])
    assert (v.diff(t)-div_stress-e["force"]).applyfunc(sp.simplify) == sp.zeros(2, 1)


def test_symbolic_both_traction_components():
    e = expressions()
    n = sp.Matrix([0, 1])
    residual = e["stress"].subs(e["z"], 0)*n+e["g"]*e["eta"]*n-e["traction"]
    assert residual.applyfunc(sp.simplify) == sp.zeros(2, 1)
    assert sp.simplify(e["traction"][0]) != 0
    assert sp.simplify(e["traction"][1]) != 0
    assert sp.simplify(e["pressure"].diff(e["x"])) != 0
