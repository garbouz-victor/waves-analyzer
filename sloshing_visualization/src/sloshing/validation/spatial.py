"""MMS errors integrated by FEM quadrature against analytic expressions."""

import numpy as np
from skfem import Basis

from ..config import SimulationConfig
from ..fem_spaces import FEMSystem
from ..time_integrator import MidpointIntegrator
from .manufactured import ManufacturedProblem


def errors(fem, integrator, state):
    p = integrator.problem
    vb = Basis(fem.mesh, fem.velocity.elem, intorder=12)
    qb = Basis(fem.mesh, fem.pressure.elem, quadrature=vb.quadrature)
    full = fem.expand_velocity(state.velocity)
    vh = vb.interpolate(full)
    qh = qb.interpolate(integrator.instantaneous_pressure(state))
    x, z = vb.global_coordinates()
    norm = lambda v: float(np.sqrt(np.sum(v*v*vb.dx)))
    v = p.evaluate("velocity", x, z, state.time)
    grad = p.evaluate("gradient", x, z, state.time)
    sb = fem.scalar.boundary("surface", intorder=12)
    coeff = np.zeros(fem.scalar.N)
    coeff[fem.surface_dofs] = state.eta
    eta_h = sb.interpolate(coeff)
    xs = sb.global_coordinates()[0]
    eta = p.evaluate("eta", xs, np.zeros_like(xs), state.time)
    vs = fem.velocity.boundary("surface", intorder=12).interpolate(full)
    qs = fem.pressure.boundary("surface", intorder=12).interpolate(integrator.instantaneous_pressure(state))
    traction = np.stack((fem.config.nu*(vs.grad[0, 1]+vs.grad[1, 0]),
                         -qs+2*fem.config.nu*vs.grad[1, 1]+fem.config.g*eta_h))
    traction -= p.surface_traction(xs, state.time)
    return {"velocity_l2": norm(vh-v), "velocity_h1": norm(vh.grad-grad),
            "pressure_l2": norm(qh-p.evaluate("pressure", x, z, state.time)),
            "eta_l2": float(np.sqrt(np.sum((eta_h-eta)**2*sb.dx))),
            "divergence_l2": norm(vh.grad[0, 0]+vh.grad[1, 1]),
            "traction_tangent_l2": float(np.sqrt(np.sum(traction[0]**2*sb.dx))),
            "traction_normal_l2": float(np.sqrt(np.sum(traction[1]**2*sb.dx)))}


def manufactured_run(n, dt=None, t_end=.1, integrator_class=MidpointIntegrator):
    dt = .0025*(8/n)**2 if dt is None else dt
    config = SimulationConfig(a=1, d=1, nx=n, nz=n, x_grading=0, z_grading=0,
                              nu=.1, t_end=t_end, dt=dt, snapshot_dt=t_end)
    fem = FEMSystem(config)
    integrator = integrator_class(fem, ManufacturedProblem(config))
    state = integrator.initial_state()
    for _ in range(config.nsteps):
        state = integrator.advance(state)
    result = {"n": n, "h": 2/n, "dt": dt, "t_end": t_end, **errors(fem, integrator, state)}
    return result


def convergence_rows(ns=(8, 16, 32), dt_factor=1.):
    rows = []
    for n in ns:
        row = manufactured_run(n, dt_factor*.0025*(8/n)**2)
        if rows:
            for key in list(row):
                if key.endswith(("_l2", "_h1")):
                    row[key+"_order"] = float(np.log(rows[-1][key]/row[key])/np.log(n/rows[-1]["n"]))
        rows.append(row)
        print(row, flush=True)
    return rows
