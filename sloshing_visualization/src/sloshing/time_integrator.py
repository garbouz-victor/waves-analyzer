"""Monolithic implicit midpoint; surface unknowns are eliminated exactly."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import splu


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


class MidpointIntegrator:
    def __init__(self, fem):
        self.fem = fem
        dt, g = fem.config.dt, fem.config.g
        lhs = fem.M / dt + fem.K / 2 + dt * g * fem.H / 4
        self.rhs_matrix = fem.M / dt - fem.K / 2 - dt * g * fem.H / 4
        # The open traction boundary fixes the pressure constant. Do NOT pin a
        # pressure dof or remove the constant test: that would lose mass balance.
        self.lu = splu(bmat([[lhs, -fem.B.T], [-fem.B, None]], format="csc"))
        self._pressure_lu = None

    def initial_state(self, eta=None):
        eta = self.fem.initial_surface() if eta is None else np.asarray(eta, dtype=float)
        if eta.shape != self.fem.surface_x.shape or not np.all(np.isfinite(eta)):
            raise ValueError("initial eta must be finite and match the surface P2 dofs")
        return State(0, 0.0, np.zeros(len(self.fem.free)), eta.copy())

    def advance(self, state):
        f = self.fem
        dt = f.config.dt
        rhs_u = self.rhs_matrix @ state.velocity - f.config.g * (f.C @ state.eta)
        rhs = np.concatenate([rhs_u, np.zeros(f.pressure.N)])
        sol = self.lu.solve(rhs)
        if not np.all(np.isfinite(sol)):
            raise FloatingPointError("Non-finite midpoint saddle-point solution")
        velocity = sol[:len(f.free)]
        vmid = (velocity + state.velocity) / 2
        eta = state.eta + dt * (f.R @ vmid)
        loss = dt * float(vmid @ (f.K @ vmid))
        previous_energy = sum(f.energy(state.velocity, state.eta))
        energy = sum(f.energy(velocity, eta))
        residual = energy - previous_energy + loss
        return State(state.step + 1, (state.step + 1) * dt, velocity, eta,
                     state.dissipated_energy + loss, residual,
                     max(state.max_step_energy_residual, abs(residual)),
                     max(state.max_step_energy_increase, energy - previous_energy))

    def instantaneous_pressure(self, state):
        """Recover q at snapshot time, including t=0 (not midpoint pressure)."""
        f = self.fem
        if self._pressure_lu is None:
            self._pressure_lu = splu(bmat([[f.M, -f.B.T], [-f.B, None]], format="csc"))
        force = -f.K @ state.velocity - f.config.g * (f.C @ state.eta)
        result = self._pressure_lu.solve(np.concatenate([force, np.zeros(f.pressure.N)]))
        if not np.all(np.isfinite(result)):
            raise FloatingPointError("Non-finite pressure recovery")
        return result[len(f.free):]
