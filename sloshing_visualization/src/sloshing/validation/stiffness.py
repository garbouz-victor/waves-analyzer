"""Scalar stability functions and fast decaying modes of the constrained FEM."""

from dataclasses import replace

import numpy as np
from scipy.linalg import eig
from scipy.sparse import bmat
from scipy.sparse.linalg import splu

from ..config import SimulationConfig
from ..fem_spaces import FEMSystem
from ..time_integrator import GAMMA, INTEGRATORS, State


def amplification(method, z):
    if method == "midpoint":
        return (1+z/2)/(1-z/2)
    return (1+(1-2*GAMMA)*z)/(1-GAMMA*z)**2


def scalar_step(method, y, lam, dt):
    if method == "midpoint":
        mid = y/(1-lam*dt/2)
        return y+dt*lam*mid
    stage1 = y/(1-GAMMA*dt*lam)
    stage2 = (y+(1-GAMMA)*dt*lam*stage1)/(1-GAMMA*dt*lam)
    return stage2


def scalar_rows():
    return [{"lambda_dt": z, "midpoint": float(amplification("midpoint", z)),
             "sdirk2": float(amplification("sdirk2", z)), "exact": float(np.exp(z))}
            for z in (-.01, -1., -100., -10000.)]


def scalar_time_rows():
    rows = []
    for method in INTEGRATORS:
        previous = None
        for dt in (.1, .05, .025, .0125):
            y = 1.
            for _ in range(round(1/dt)):
                y = scalar_step(method, y, -1., dt)
            error = abs(y-np.exp(-1))
            row = {"integrator": method, "dt": dt, "error": error}
            if previous is not None:
                row["order"] = float(np.log2(previous/error))
            rows.append(row)
            previous = error
    return rows


def semidiscrete_generator(fem):
    f = fem
    nv, ns = len(f.free), len(f.surface_x)
    lu = splu(bmat([[f.M, -f.B.T], [-f.B, None]], format="csc"))
    rhs = np.hstack((-f.K.toarray(), -f.config.g*f.C.toarray()))
    acceleration = lu.solve(np.vstack((rhs, np.zeros((f.pressure.N, nv+ns)))))[:nv]
    return np.vstack((acceleration, np.hstack((f.R.toarray(), np.zeros((ns, ns))))))


def stiff_fem_comparison():
    config = SimulationConfig(nx=4, nz=6, nu=1., t_end=0)
    fem = FEMSystem(config)
    values, vectors = eig(semidiscrete_generator(fem))
    candidates = np.flatnonzero((values.real < -1) & (abs(values.imag) < 1e-7*abs(values.real)))
    i = candidates[np.argmin(values[candidates].real)]
    lam = float(values[i].real)
    vec = vectors[:, i].real
    nv = len(fem.free)
    vec /= np.sqrt(2*sum(fem.energy(vec[:nv], vec[nv:])))
    dt = 100/abs(lam)
    # The FEM geometry and spatial matrices do not depend on dt.
    fem.config = replace(config, dt=dt, snapshot_dt=dt, t_end=5*dt)
    rows = []
    for method, cls in INTEGRATORS.items():
        it = cls(fem)
        state = State(0, 0, vec[:nv].copy(), vec[nv:].copy())
        for step in range(1, 6):
            state = it.advance(state)
            exact_factor = float(np.exp(lam*state.time))
            error = np.sqrt(2*sum(fem.energy(state.velocity-exact_factor*vec[:nv],
                                             state.eta-exact_factor*vec[nv:])))
            amplitude = float(vec[:nv]@(fem.M@state.velocity)+fem.config.g*vec[nv:]@(fem.S@state.eta))
            rows.append({"integrator": method, "lambda": lam, "dt": dt, "lambda_dt": lam*dt,
                         "step": step, "time": state.time, "exact_decay": exact_factor,
                         "amplitude": amplitude, "energy_norm_error": float(error),
                         "weak_constraint_max": float(np.max(abs(fem.B@state.velocity)))})
    return rows
