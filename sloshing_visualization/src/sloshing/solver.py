"""Solver facade: no plotting, particles, or fabricated mode damping."""

from .diagnostics import Diagnostics
from .fem_spaces import FEMSystem
from .time_integrator import INTEGRATORS


class SloshingSolver:
    def __init__(self, config):
        self.config = config
        self.fem = FEMSystem(config)
        self.integrator = INTEGRATORS[config.integrator](self.fem)

    def snapshots(self, eta_initial=None):
        state = self.integrator.initial_state(eta_initial)
        diagnostics = Diagnostics(self.fem, state)
        row = diagnostics.measure(state)
        diagnostics.validate(row)
        yield state, row
        for _ in range(self.config.nsteps):
            state = self.integrator.advance(state)
            if state.step % self.config.save_every == 0 or state.step == self.config.nsteps:
                row = diagnostics.measure(state)
                diagnostics.validate(row)
                yield state, row
