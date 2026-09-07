"""Impermeability is essential; tangential Navier traction is natural."""
import numpy as np
from dolfinx import fem
from petsc4py import PETSc


def velocity_pressure_bcs(space, config, facets):
    bcs = []
    fdim = space.mesh.topology.dim-1
    for side, indices in facets.items():
        components = (0, 1) if side in config.no_slip_walls else (
            (0,) if side in ("left", "right") else (1,))
        for component in components:
            subspace = space.sub(0).sub(component)
            dofs = fem.locate_dofs_topological(subspace, fdim, indices)
            bcs.append(fem.dirichletbc(PETSc.ScalarType(0), dofs, subspace))
    # One pressure DOF removes the redundant constant, not a compressibility penalty.
    pspace, _ = space.sub(1).collapse()
    pdofs = fem.locate_dofs_geometrical((space.sub(1), pspace), lambda x:
        np.isclose(x[0], config.x_min) & np.isclose(x[1], config.z_min))
    zero = fem.Function(pspace)
    bcs.append(fem.dirichletbc(zero, pdofs, space.sub(1)))
    return bcs
