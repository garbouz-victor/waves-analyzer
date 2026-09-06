import numpy as np

from sloshing.config import SimulationConfig
from sloshing.solver import SloshingSolver
from sloshing.validation.qualification.fem_evaluation import FEMGeometry


def test_fem_particle_evaluation_preserves_wall_zero_and_limit():
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, alpha_deg=.02, integrator="sdirk2", t_end=.05))
    f, it = solver.fem, solver.integrator
    state = it.initial_state()
    for _ in range(f.config.nsteps):
        state = it.advance(state)
    full = f.expand_velocity(state.velocity)
    geom = FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs)
    for sign in (-1, 1):
        eps = np.array([.01, .001, .0001, .00001, 0.])
        for z in (-.05, -.2, -1.):
            points = np.vstack((sign*(1-eps), np.full(len(eps), z)))
            v = np.stack([geom.evaluate(full[ids], points) for ids in f.component_dofs])
            speed = np.linalg.norm(v, axis=0)
            assert speed[-1] < 1e-15
            assert speed[-2] < .01*speed[0]
