"""Declared finite-slip linear FEM. Historical no-slip classes stay unchanged."""
from dataclasses import dataclass
import numpy as np
from skfem import BilinearForm, asm

from ..config import SimulationConfig
from ..diagnostics import Diagnostics
from ..fem_spaces import FEMSystem
from ..time_integrator import SDIRK2Integrator, State, GAMMA

LABEL = "DECLARED_EXTENSION — конечное скольжение Навье"


@dataclass(frozen=True)
class NavierConfig(SimulationConfig):
    slip_length_m: float = .50

    def __post_init__(self):
        super().__post_init__()
        if self.slip_length_m <= 0:
            raise ValueError("Finite positive side-wall slip length required")
        if self.integrator != "sdirk2":
            raise ValueError("Declared extension currently uses SDIRK2")


@BilinearForm
def tangential_wall_mass(u, v, _):
    return u[1] * v[1]


class NavierFEM(FEMSystem):
    def __init__(self, config):
        super().__init__(config)
        horizontal = self.component_dofs[0]
        normals = np.intersect1d(np.r_[self.wall_dofs["left"], self.wall_dofs["right"]], horizontal)
        self.fixed = np.unique(np.r_[normals, self.wall_dofs["bottom"]])
        self.free = np.setdiff1d(np.arange(self.velocity.N), self.fixed)
        self.side_basis = {name: self.velocity.boundary(name, intorder=6) for name in ("left", "right")}
        self.K_bulk_full = self.K_full.copy()
        self.K_wall_full = sum(asm(tangential_wall_mass, basis) for basis in self.side_basis.values()).tocsr() * (config.nu/config.slip_length_m)
        self.K_full = self.K_bulk_full + self.K_wall_full
        self.M = self.M_full[self.free][:, self.free].tocsr()
        self.K_bulk = self.K_bulk_full[self.free][:, self.free].tocsr()
        self.K_wall = self.K_wall_full[self.free][:, self.free].tocsr()
        self.K = self.K_bulk + self.K_wall
        self.B = self.B_full[:, self.free].tocsr()
        self.R = self.R_full[:, self.free].tocsr()
        self.C = (self.R.T @ self.S).tocsr()
        self.H = (self.C @ self.R).tocsr()
        if [self.R.getrow(i).nnz for i in (0, self.R.shape[0]-1)] != [1, 1]:
            raise ValueError("Vertical endpoint is not free")
        self.surface_basis = self.velocity.boundary("surface", intorder=6)


@dataclass
class NavierState(State):
    wall_dissipated_energy: float = 0.
    bulk_dissipated_energy: float = 0.


class NavierIntegrator(SDIRK2Integrator):
    def initial_state(self, eta=None):
        return NavierState(**super().initial_state(eta).__dict__)

    def _stage(self, *args):
        result = super()._stage(*args)
        self.stages.append(result)
        return result

    def advance(self, state):
        self.stages = []
        result = super().advance(state)
        f, dt = self.fem, self.fem.config.dt
        weights = (1-GAMMA, GAMMA)
        wall = dt*sum(b*float(s[0] @ (f.K_wall @ s[0])) for b, s in zip(weights, self.stages))
        bulk = dt*sum(b*float(s[0] @ (f.K_bulk @ s[0])) for b, s in zip(weights, self.stages))
        self.wall_step_loss, self.bulk_step_loss = wall, bulk
        residuals = []
        for v, eta, kv, _, load, p in self.stages:
            terms = (f.M @ kv, f.K @ v, f.config.g*(f.C @ eta), -(f.B.T @ p), -load)
            residuals.append(float(np.linalg.norm(sum(terms))) / max(1e-30, sum(np.linalg.norm(t) for t in terms)))
        self.momentum_residual = max(residuals)
        return NavierState(**result.__dict__,
                           wall_dissipated_energy=state.wall_dissipated_energy+wall,
                           bulk_dissipated_energy=state.bulk_dissipated_energy+bulk)


class NavierDiagnostics(Diagnostics):
    def measure(self, state):
        row = super().measure(state)
        f, c = self.fem, self.fem.config
        full = f.expand_velocity(state.velocity)
        row.update(max_eta_over_a=float(np.max(abs(state.eta)))/c.a,
                   max_eta_over_b=float(np.max(abs(state.eta)))/c.slip_length_m,
                   wall_dissipation=float(state.velocity @ (f.K_wall @ state.velocity)),
                   bulk_dissipation=float(state.velocity @ (f.K_bulk @ state.velocity)),
                   cumulative_wall_dissipation=state.wall_dissipated_energy,
                   cumulative_bulk_dissipation=state.bulk_dissipated_energy)
        physical_scale = c.g*max(abs(c.slope), 1e-12)
        Uscale = np.sqrt(c.g*c.a)*max(abs(c.slope), 1e-12)
        field = f.velocity.interpolate(full)
        adv = np.einsum("ij...,j...->i...", field.grad, np.asarray(field))
        row["convective_indicator"] = float(np.max(np.sqrt(np.sum(adv**2, axis=0)))) / physical_scale
        trace = f.surface_basis.interpolate(full)
        scalar_trace = f.scalar.boundary("surface", intorder=6)
        eta_scalar = np.zeros(f.scalar.N)
        eta_scalar[f.surface_dofs] = state.eta
        eta_field = scalar_trace.interpolate(eta_scalar)
        omitted = np.asarray(eta_field)*trace.grad[1, 1] - np.asarray(trace)[0]*eta_field.grad[0]
        row["kinematic_indicator"] = float(np.max(abs(omitted))) / Uscale
        wall_power, wall_residual = 0., 0.
        for name, normal in (("left", -1.), ("right", 1.)):
            basis = f.side_basis[name]
            value = basis.interpolate(full)
            traction = c.nu*(normal*(value.grad[1, 0]+value.grad[0, 1]) + np.asarray(value)[1]/c.slip_length_m)
            wall_residual += float(np.sum(traction**2 * basis.dx))
            wall_power += c.nu/c.slip_length_m*float(np.sum(np.asarray(value)[1]**2*basis.dx))
        row["wall_traction_L2"] = np.sqrt(wall_residual)
        row["wall_power_quadrature_error"] = abs(wall_power-row["wall_dissipation"])
        return row

    def validate(self, row):
        if not all(np.isfinite(v) for v in row.values()):
            raise ValueError("Nonfinite Navier fields/diagnostics")
        limits = {"weak_divergence_l2": 1e-9, "volume_error": 1e-10,
                  "energy_balance_relative": 1e-8, "left_u_max": 1e-12, "right_u_max": 1e-12,
                  "bottom_u_max": 1e-12, "bottom_w_max": 1e-12, "wall_power_quadrature_error": 1e-11}
        bad = {k: row[k] for k, tolerance in limits.items() if abs(row[k]) > tolerance}
        if bad or min(row["wall_dissipation"], row["bulk_dissipation"]) < -1e-15:
            raise ValueError("Navier physical budget/BC violation: " + repr(bad))


def channel_benchmark(config):
    """Steady SIDE-WALL channel only; no vessel bottom/free surface in this test."""
    from scipy.sparse.linalg import spsolve
    f = NavierFEM(config)
    wdofs = f.component_dofs[1]
    G = .001
    lhs = (f.K_bulk_full + f.K_wall_full)[wdofs][:, wdofs]
    rhs = G*(f.M_full[wdofs][:, wdofs] @ np.ones(len(wdofs)))
    w = spsolve(lhs, rhs)
    x = f.velocity.doflocs[0, wdofs]
    exact = G/(2*config.nu)*(config.a**2-x*x+2*config.a*config.slip_length_m)
    return {"benchmark": "steady_vertical_channel_with_Navier_sides_only",
            "max_absolute_error_m_s": float(np.max(abs(w-exact))),
            "relative_error": float(np.max(abs(w-exact))/np.max(abs(exact))),
            "exact_formula": "G/(2*nu)*(a*a-x*x+2*a*b)", "G_m_s2": G}
