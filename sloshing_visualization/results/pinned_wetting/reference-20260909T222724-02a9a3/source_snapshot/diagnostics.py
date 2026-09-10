"""Physical norms and exact discrete midpoint/SDIRK2 energy budgets."""

import warnings

import numpy as np


PINNED_MESSAGE = "Pinned contact line: consequence of strict no-slip model"
SLOPE_WARNING = "Local small-slope assumption is breaking down near the wall."


class PhysicsViolation(RuntimeError):
    pass


def surface_slopes(x, eta):
    """Exact one-sided extrema of each P2 edge derivative, no finite differences."""
    left, middle, right = eta[:-2:2], eta[1::2], eta[2::2]
    h = x[2::2] - x[:-2:2]
    return (-3 * left + 4 * middle - right) / h, (left - 4 * middle + 3 * right) / h


def quadratic_edge_max(left, middle, right):
    """Exact max |P2| on boundary edges, including interior extrema."""
    aa, bb = 2*left - 4*middle + 2*right, -3*left + 4*middle - right
    stationary = np.divide(-bb, 2*aa, out=np.zeros_like(aa), where=aa != 0)
    interior = (aa != 0) & (stationary > 0) & (stationary < 1)
    values = np.where(interior, aa*stationary**2 + bb*stationary + left, left)
    return float(max(np.max(abs(left)), np.max(abs(right)), np.max(abs(values))))


class Diagnostics:
    def __init__(self, fem, initial):
        self.fem = fem
        self.initial_energy = sum(fem.energy(initial.velocity, initial.eta))
        self.initial_volume = float(fem.surface_weights @ initial.eta)
        self.contacts = initial.eta[[0, -1]].copy()
        self.slope_warned = False
        self.divergence_warned = False
        self.wall_scalar_dofs = {}
        for name in fem.wall_dofs:
            facets = fem.mesh.boundaries[name]
            vertices = fem.mesh.facets[:, facets]
            self.wall_scalar_dofs[name] = (fem.scalar.nodal_dofs[0, vertices[0]],
                                          fem.scalar.facet_dofs[0, facets],
                                          fem.scalar.nodal_dofs[0, vertices[1]])

    def measure(self, state):
        f = self.fem
        velocity = f.expand_velocity(state.velocity)
        field = f.velocity.interpolate(velocity)
        gradient = field.grad
        divergence = gradient[0, 0] + gradient[1, 1]
        div_l2 = float(np.sqrt(np.sum(divergence**2 * f.velocity.dx)))
        grad_l2 = float(np.sqrt(np.sum(gradient**2 * f.velocity.dx[None, None, :, :])))
        weak = f.B @ state.velocity
        weak_l2 = float(np.sqrt(max(0., weak @ f.pressure_mass_lu.solve(weak))))
        kinetic, potential = f.energy(state.velocity, state.eta)
        slopes_left, slopes_right = surface_slopes(f.surface_x, state.eta)
        max_slope = float(max(np.max(abs(slopes_left)), np.max(abs(slopes_right))))
        volume = float(f.surface_weights @ state.eta)
        energy_scale = max(self.initial_energy, 1e-30)
        budget = (kinetic+potential+state.dissipated_energy+state.rk_energy_correction
                  -state.work_energy-self.initial_energy)
        result = {
            "time": state.time,
            "kinetic_energy": kinetic,
            "gravitational_potential_energy": potential,
            "total_energy": kinetic + potential,
            "viscous_dissipation": float(state.velocity @ (f.K @ state.velocity)),
            "cumulative_viscous_dissipation": state.dissipated_energy,
            "cumulative_work": state.work_energy,
            "rk_energy_correction": state.rk_energy_correction,
            "energy_balance_residual": budget,
            "energy_balance_relative": budget / energy_scale,
            "step_energy_residual": state.step_energy_residual,
            "max_step_energy_residual": state.max_step_energy_residual,
            "max_step_energy_increase": state.max_step_energy_increase,
            "divergence_l2": div_l2,
            "velocity_gradient_l2": grad_l2,
            "divergence_relative": div_l2 / max(grad_l2, 1e-30),
            "weak_divergence_l2": weak_l2,
            "surface_integral": volume,
            "volume_error": abs(volume - self.initial_volume),
            "eta_left": float(state.eta[0]),
            "eta_right": float(state.eta[-1]),
            "contact_error": float(np.max(abs(state.eta[[0, -1]] - self.contacts))),
            "max_slope": max_slope,
            "slope_left": float(slopes_left[0]),
            "slope_right": float(slopes_right[-1]),
        }
        for wall, (left, middle, right) in self.wall_scalar_dofs.items():
            for component, comp_dofs in zip(("u", "w"), f.component_dofs):
                values = velocity[comp_dofs]
                result[f"{wall}_{component}_max"] = quadratic_edge_max(values[left], values[middle], values[right])
        return result

    def validate(self, row):
        c = self.fem.config
        if not all(np.isfinite(value) for value in row.values()):
            raise PhysicsViolation(f"Non-finite diagnostics at t={row['time']}")
        checks = {
            "weak_divergence_l2": c.weak_divergence_tolerance,
            "volume_error": c.volume_tolerance,
            "contact_error": c.contact_tolerance,
            "energy_balance_relative": c.energy_relative_tolerance,
            "max_step_energy_residual": c.energy_relative_tolerance * max(self.initial_energy, 1e-30),
            "max_step_energy_increase": c.energy_relative_tolerance * max(self.initial_energy, 1e-30),
        }
        checks.update({f"{wall}_{component}_max": c.wall_tolerance
                       for wall in ("left", "right", "bottom") for component in ("u", "w")})
        violations = [f"{key}={row[key]:.4e} > {limit:.4e}"
                      for key, limit in checks.items() if abs(row[key]) > limit]
        if violations:
            raise PhysicsViolation(f"t={row['time']}: " + "; ".join(violations))
        if row["max_slope"] > c.small_slope_limit and not self.slope_warned:
            warnings.warn(SLOPE_WARNING + f" max |eta_x|={row['max_slope']:.4g}", RuntimeWarning, stacklevel=2)
            self.slope_warned = True
        if (row["divergence_relative"] > c.divergence_relative_warning
                and row["divergence_l2"] > c.weak_divergence_tolerance
                and not self.divergence_warned):
            warnings.warn("Strong divergence is under-resolved: "
                          f"||div v||/||grad v||={row['divergence_relative']:.3g}; refine the mesh.",
                          RuntimeWarning, stacklevel=2)
            self.divergence_warned = True
