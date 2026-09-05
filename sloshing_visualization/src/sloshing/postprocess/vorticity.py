"""Exact elementwise P1 curl of P2 velocity, without cross-element smoothing."""

from scipy.sparse.linalg import splu
from skfem import Basis, BilinearForm, ElementTriP1DG, asm

from ..fem_spaces import scalar_mass


@BilinearForm
def curl_form(u, test, _):
    return (u.grad[1, 0] - u.grad[0, 1]) * test


class VorticityProjector:
    def __init__(self, fem):
        self.basis = Basis(fem.mesh, ElementTriP1DG(), quadrature=fem.velocity.quadrature)
        self.mass_lu = splu(asm(scalar_mass, self.basis).tocsc())
        self.curl = asm(curl_form, fem.velocity, self.basis).tocsr()

    def evaluate(self, full_velocity):
        return self.mass_lu.solve(self.curl @ full_velocity)
