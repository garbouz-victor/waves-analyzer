"""Monolithic second-order methods for the constrained velocity/surface DAE."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import splu

from .problems import ProblemLoads, ProductionProblem


GAMMA = 1 - 1/np.sqrt(2.)


@dataclass
class State:
    step: int
    time: float
    velocity: np.ndarray
    eta: np.ndarray
    dissipated_energy: float = 0.0
    step_energy_residual: float = 0.0
    max_step_energy_residual: float = 0.0
    max_step_energy_increase: float = 0.0
    work_energy: float = 0.0
    rk_energy_correction: float = 0.0


class IntegratorBase:
    def __init__(self, fem, problem=None):
        self.fem = fem
        self.problem = ProductionProblem(fem.config) if problem is None else problem
        self.loads = ProblemLoads(fem, self.problem)
        self._pressure_lu = None

    def initial_state(self, eta=None):
        velocity, default_eta = self.loads.initial_fields()
        eta = default_eta if eta is None else np.asarray(eta, dtype=float)
        if eta.shape != self.fem.surface_x.shape or not np.all(np.isfinite(eta)):
            raise ValueError("initial eta must be finite and match the surface P2 dofs")
        return State(0, 0.0, velocity, eta.copy())

    def _finish(self, state, velocity, eta, loss, work=0., correction=0.):
        energy = sum(self.fem.energy(velocity, eta))
        previous = sum(self.fem.energy(state.velocity, state.eta))
        residual = energy-previous+loss-work+correction
        return State(state.step+1, (state.step+1)*self.fem.config.dt, velocity, eta,
                     state.dissipated_energy+loss, residual,
                     max(state.max_step_energy_residual, abs(residual)),
                     max(state.max_step_energy_increase, energy-previous),
                     state.work_energy+work, state.rk_energy_correction+correction)

    def instantaneous_pressure(self, state):
        """q at snapshot time (not stage pressure), with its traction-fixed gauge."""
        f = self.fem
        if self._pressure_lu is None:
            self._pressure_lu = splu(bmat([[f.M, -f.B.T], [-f.B, None]], format="csc"))
        force = -f.K@state.velocity-f.config.g*(f.C@state.eta)+self.loads(state.time)
        result = self._pressure_lu.solve(np.r_[force, np.zeros(f.pressure.N)])
        if not np.all(np.isfinite(result)):
            raise FloatingPointError("Non-finite pressure recovery")
        return result[len(f.free):]


class MidpointIntegrator(IntegratorBase):
    def __init__(self, fem, problem=None):
        super().__init__(fem, problem)
        dt, g = fem.config.dt, fem.config.g
        lhs = fem.M/dt+fem.K/2+dt*g*fem.H/4
        self.rhs_matrix = fem.M/dt-fem.K/2-dt*g*fem.H/4
        self.lu = splu(bmat([[lhs, -fem.B.T], [-fem.B, None]], format="csc"))

    def advance(self, state):
        f, dt = self.fem, self.fem.config.dt
        load = self.loads(state.time+dt/2)
        rhs = self.rhs_matrix@state.velocity-f.config.g*(f.C@state.eta)+load
        sol = self.lu.solve(np.r_[rhs, np.zeros(f.pressure.N)])
        if not np.all(np.isfinite(sol)):
            raise FloatingPointError("Non-finite midpoint saddle-point solution")
        velocity = sol[:len(f.free)]
        vmid = (velocity+state.velocity)/2
        eta = state.eta+dt*(f.R@vmid)
        return self._finish(state, velocity, eta, dt*float(vmid@(f.K@vmid)),
                            dt*float(vmid@load))


class SDIRK2Integrator(IntegratorBase):
    """Stiffly accurate SDIRK2, gamma=1-1/sqrt(2).

    Each stage solves [M/s+K+s*g*R.T*S*R, -B.T; -B, 0], s=gamma*dt.
    Stage eta=eta_base+s*R*V is an exact elimination, not an explicit update.
    """

    def __init__(self, fem, problem=None):
        super().__init__(fem, problem)
        s = GAMMA*fem.config.dt
        lhs = fem.M/s+fem.K+s*fem.config.g*fem.H
        self.lu = splu(bmat([[lhs, -fem.B.T], [-fem.B, None]], format="csc"))

    def _stage(self, vbase, ebase, t):
        f, s = self.fem, GAMMA*self.fem.config.dt
        load = self.loads(t)
        rhs = f.M@vbase/s-f.config.g*(f.C@ebase)+load
        sol = self.lu.solve(np.r_[rhs, np.zeros(f.pressure.N)])
        if not np.all(np.isfinite(sol)):
            raise FloatingPointError("Non-finite SDIRK2 stage")
        v = sol[:len(f.free)]
        eta = ebase+s*(f.R@v)
        return v, eta, (v-vbase)/s, f.R@v, load, sol[len(f.free):]

    def advance(self, state):
        f, dt = self.fem, self.fem.config.dt
        b1, b2 = 1-GAMMA, GAMMA
        v1, e1, kv1, ke1, load1, _ = self._stage(state.velocity, state.eta, state.time+GAMMA*dt)
        v2, e2, kv2, ke2, load2, pressure2 = self._stage(state.velocity+dt*b1*kv1,
                                              state.eta+dt*b1*ke1, state.time+dt)
        # Stiff accuracy: stage 2 IS the endpoint state. B*k_v2=0, hence
        # its pressure also solves the instantaneous acceleration saddle system.
        # Exposing it avoids a redundant LU solve for each saved snapshot.
        self.last_stage_pressure = pressure2
        loss = dt*(b1*float(v1@(f.K@v1))+b2*float(v2@(f.K@v2)))
        work = dt*(b1*float(v1@load1)+b2*float(v2@load2))
        # RK identity: m_ij=b_i*a_ij+b_j*a_ji-b_i*b_j=diag(-gamma²,+gamma²).
        # The correction is SIGNED: SDIRK2 is not algebraically stable, and
        # this term must not be described as physical heat.
        correction = dt**2*GAMMA**2*(sum(f.energy(kv2, ke2))-sum(f.energy(kv1, ke1)))
        return self._finish(state, v2, e2, loss, work, correction)


INTEGRATORS = {"midpoint": MidpointIntegrator, "sdirk2": SDIRK2Integrator}
