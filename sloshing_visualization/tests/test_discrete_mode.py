"""A viscous FEM eigenmode, not an inviscid potential-flow frequency."""

import numpy as np
from scipy.linalg import eig
from scipy.sparse import bmat
from scipy.sparse.linalg import splu

from sloshing import SimulationConfig, SloshingSolver
from sloshing.time_integrator import State


def test_midpoint_against_independent_viscous_semidiscrete_eigenmode():
    solver = SloshingSolver(SimulationConfig(nx=4, nz=4, dt=.001, snapshot_dt=.01, t_end=.05))
    f = solver.fem
    nv, ns = len(f.free), len(f.surface_x)
    mass_stokes = splu(bmat([[f.M, -f.B.T], [-f.B, None]], format="csc"))
    forces = np.hstack((-f.K.toarray(), -f.config.g*f.C.toarray()))
    accelerations = mass_stokes.solve(np.vstack((forces, np.zeros((f.pressure.N, nv+ns)))))[:nv]
    generator = np.vstack((accelerations, np.hstack((f.R.toarray(), np.zeros((ns, ns))))))
    eigenvalues, eigenvectors = eig(generator)
    candidates = np.flatnonzero((eigenvalues.imag > .1) & (eigenvalues.real < -1e-5))
    assert len(candidates) > 0
    mode = candidates[np.argmin(abs(eigenvalues[candidates].imag - 4.0))]
    lam, vec = eigenvalues[mode], eigenvectors[:, mode]
    np.testing.assert_allclose(f.B @ vec[:nv], 0, atol=1e-10)
    np.testing.assert_allclose(vec[nv:][[0, -1]], 0, atol=1e-10)
    state = State(0, 0, vec[:nv].real.copy(), vec[nv:].real.copy())
    e0 = sum(f.energy(state.velocity, state.eta))
    for _ in range(solver.config.nsteps):
        state = solver.integrator.advance(state)
    exact = (np.exp(lam * state.time) * vec).real
    error = sum(f.energy(state.velocity-exact[:nv], state.eta-exact[nv:]))
    assert np.sqrt(error/e0) < 5e-6
    assert sum(f.energy(state.velocity, state.eta)) < e0
