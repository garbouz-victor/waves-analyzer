"""P2/P1 Taylor–Hood and the exact P2 velocity trace on the free surface."""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from skfem import Basis, BilinearForm, ElementTriP1, ElementTriP2, ElementVector, asm
from skfem.helpers import ddot, div, dot, sym_grad

from .mesh import generate_mesh


@BilinearForm
def velocity_mass(u, v, _):
    return dot(u, v)


@BilinearForm
def strain_form(u, v, _):
    return 2.0 * ddot(sym_grad(u), sym_grad(v))


@BilinearForm
def divergence_form(u, p, _):
    return p * div(u)


@BilinearForm
def scalar_mass(u, v, _):
    return u * v


class FEMSystem:
    def __init__(self, config):
        self.config = config
        self.mesh = generate_mesh(config)
        self.velocity = Basis(self.mesh, ElementVector(ElementTriP2()), intorder=6)
        self.scalar = Basis(self.mesh, ElementTriP2(), intorder=6)
        self.pressure = Basis(self.mesh, ElementTriP1(), intorder=6)
        self.component_dofs = self.velocity.split_indices()
        self.wall_dofs = {name: self.velocity.get_dofs(name).all()
                          for name in ("left", "right", "bottom")}
        self.fixed = np.unique(np.concatenate(list(self.wall_dofs.values())))
        self.free = np.setdiff1d(np.arange(self.velocity.N), self.fixed)
        self.M_full = asm(velocity_mass, self.velocity).tocsr()
        self.K_full = config.nu * asm(strain_form, self.velocity).tocsr()
        self.B_full = asm(divergence_form, self.velocity, self.pressure).tocsr()
        self.M = self.M_full[self.free][:, self.free].tocsr()
        self.K = self.K_full[self.free][:, self.free].tocsr()
        self.B = self.B_full[:, self.free].tocsr()
        self.pressure_mass = asm(scalar_mass, self.pressure).tocsc()
        self.pressure_mass_lu = splu(self.pressure_mass)
        top = self.scalar.get_dofs("surface").all()
        self.surface_dofs = top[np.argsort(self.scalar.doflocs[0, top])]
        self.surface_x = self.scalar.doflocs[0, self.surface_dofs]
        sbasis = self.scalar.boundary("surface", intorder=6)
        surface_mass = asm(scalar_mass, sbasis).tocsr()
        self.S = surface_mass[self.surface_dofs][:, self.surface_dofs].tocsr()
        # R includes both endpoints; their velocity columns are eliminated.
        ns = len(top)
        self.R_full = csr_matrix((np.ones(ns), (np.arange(ns),
                                  self.component_dofs[1][self.surface_dofs])),
                                 shape=(ns, self.velocity.N))
        self.R = self.R_full[:, self.free].tocsr()
        self.C = (self.R.T @ self.S).tocsr()
        self.H = (self.C @ self.R).tocsr()
        self.surface_weights = np.asarray(self.S @ np.ones(ns)).ravel()

    def expand_velocity(self, reduced):
        full = np.zeros(self.velocity.N)
        full[self.free] = reduced
        return full

    def initial_surface(self):
        return self.config.slope * self.surface_x

    def energy(self, velocity, eta):
        kinetic = 0.5 * float(velocity @ (self.M @ velocity))
        potential = 0.5 * self.config.g * float(eta @ (self.S @ eta))
        return kinetic, potential

    def velocity_at(self, full_velocity, points):
        probe = self.scalar.probes(np.asarray(points))
        return np.stack([probe @ full_velocity[dofs] for dofs in self.component_dofs])
