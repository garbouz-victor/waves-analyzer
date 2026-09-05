"""Analytic MMS derived with SymPy, independently of all FEM matrices."""

from functools import lru_cache

import numpy as np
import sympy as sp


@lru_cache(None)
def expressions():
    x, z, t = sp.symbols("x z t", real=True)
    nu, g = sp.symbols("nu g", positive=True)
    X, Z, A = (1-x**2)**2, (z+1)**2, sp.exp(-t)
    psi = A*X*Z
    v = sp.Matrix([sp.diff(psi, z), -sp.diff(psi, x)])
    eta, q = A*sp.diff(X, x), A*x*(z+1)
    gradient = v.jacobian([x, z])
    D = (gradient + gradient.T)/2
    stress = -q*sp.eye(2) + 2*nu*D
    laplacian = v.diff(x, 2)+v.diff(z, 2)
    force = sp.simplify(v.diff(t)+sp.Matrix([q.diff(x), q.diff(z)])-nu*laplacian)
    normal = sp.Matrix([0, 1])
    traction = sp.simplify(stress.subs(z, 0)*normal+g*eta*normal)
    return dict(x=x, z=z, t=t, nu=nu, g=g, psi=psi, velocity=v, eta=eta,
                pressure=q, gradient=gradient, stress=stress, force=force, traction=traction)


class ManufacturedProblem:
    name = "manufactured_exp_streamfunction"
    is_production = False

    @staticmethod
    def load_time_factor(t):
        """Analytic separability f(t)=exp(-t)f(0), r(t)=exp(-t)r(0)."""
        return np.exp(-t)

    def __init__(self, config):
        if config.a != 1 or config.d != 1 or config.x_grading or config.z_grading:
            raise ValueError("MMS requires a=d=1 and uniform coordinates; production geometry is unchanged")
        self.config = config
        e = expressions()
        self.functions = {}
        for key in ("velocity", "gradient", "eta", "pressure", "force", "traction"):
            expr = e[key].subs({e["nu"]: config.nu, e["g"]: config.g})
            shape = expr.shape if isinstance(expr, sp.MatrixBase) else ()
            if shape == (2, 1):
                shape = (2,)
            entries = list(expr) if shape else [expr]
            self.functions[key] = (shape, [sp.lambdify((e["x"], e["z"], e["t"]), entry, "numpy")
                                           for entry in entries])

    def evaluate(self, key, x, z, t):
        x, z = np.broadcast_arrays(x, z)
        shape, funcs = self.functions[key]
        values = [np.broadcast_to(fn(x, z, t), x.shape) for fn in funcs]
        return np.stack(values).reshape(shape+x.shape) if shape else values[0]

    def initial_velocity(self, x, z):
        return self.evaluate("velocity", x, z, 0.)

    def initial_eta(self, x):
        return self.evaluate("eta", x, np.zeros_like(x), 0.)

    def body_force(self, x, z, t):
        return self.evaluate("force", x, z, t)

    def surface_traction(self, x, t):
        return self.evaluate("traction", x, np.zeros_like(x), t)
