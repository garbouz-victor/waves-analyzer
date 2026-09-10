"""Zero-source production definition and optional analytic validation loads."""

import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import splu
from skfem import Basis, LinearForm, asm
from skfem.helpers import dot


class ProductionProblem:
    name = "physical_release"
    is_production = True

    def __init__(self, config):
        self.config = config

    def initial_velocity(self, x, z):
        return np.zeros((2,) + np.broadcast(x, z).shape)

    def initial_eta(self, x):
        return x*self.config.slope

    def body_force(self, x, z, t):
        return np.zeros((2,) + np.broadcast(x, z).shape)

    def surface_traction(self, x, t):
        return np.zeros((2,) + np.shape(x))


@LinearForm
def vector_load(v, w):
    return dot(w.load, v)


@LinearForm
def scalar_load(v, w):
    return w.load*v


class ProblemLoads:
    def __init__(self, fem, problem):
        self.fem, self.problem = fem, problem
        self.zero = np.zeros(len(fem.free))
        self.cached = None
        if not problem.is_production:
            # MMS load times P2 has total degree up to 7; degree 6 is insufficient.
            self.volume = Basis(fem.mesh, fem.velocity.elem, intorder=10)
            self.surface = self.volume.boundary("surface", intorder=10)
            if hasattr(problem, "load_time_factor"):
                self.cached = self._assemble(0.)

    def _assemble(self, t):
        p = self.problem
        x, z = self.volume.global_coordinates()
        xs = self.surface.global_coordinates()[0]
        full = asm(vector_load, self.volume, load=p.body_force(x, z, t))
        full += asm(vector_load, self.surface, load=p.surface_traction(xs, t))
        return full[self.fem.free]

    def __call__(self, t):
        if self.problem.is_production:
            return self.zero
        if self.cached is not None:
            return self.problem.load_time_factor(t)*self.cached
        return self._assemble(t)

    def initial_fields(self):
        f, p = self.fem, self.problem
        if p.is_production:
            return np.zeros(len(f.free)), f.initial_surface()
        x, z = self.volume.global_coordinates()
        rhs = asm(vector_load, self.volume, load=p.initial_velocity(x, z))[f.free]
        saddle = splu(bmat([[f.M, -f.B.T], [-f.B, None]], format="csc"))
        # Constrained L2 projection: do not initialize an inconsistent DAE state.
        velocity = saddle.solve(np.r_[rhs, np.zeros(f.pressure.N)])[:len(f.free)]
        sb = f.scalar.boundary("surface", intorder=10)
        eta_rhs = asm(scalar_load, sb, load=p.initial_eta(sb.global_coordinates()[0]))[f.surface_dofs]
        eta = np.zeros_like(f.surface_x)
        eta[[0, -1]] = p.initial_eta(f.surface_x[[0, -1]])
        eta[1:-1] = splu(f.S[1:-1, 1:-1].tocsc()).solve((eta_rhs-f.S@eta)[1:-1])
        return velocity, eta
