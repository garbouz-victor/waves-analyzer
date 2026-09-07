"""Dimensional, predeclared STEP 3A.1 policies; never tuned in benchmark logic."""
from dataclasses import dataclass


@dataclass(frozen=True)
class EnergyPolicy:
    mass_relative_tolerance: float = 1e-8
    growth_absolute_tolerance: float = 1e-9  # J/m
    budget_relative_tolerance: float = .01
    change_relative_tolerance: float = .05
    change_floor_absolute: float = 1e-10  # J/m
    change_floor_fraction: float = 1e-6
    scale_floor: float = 1e-30
    negative_power_tolerance: float = 1e-12  # W/m
    interval_fraction_warning: float = .5


@dataclass(frozen=True)
class SettlingPolicy:
    window_s: float = .10
    max_angle_rate_deg_per_s: float = .5
    max_contact_line_speed_m_per_s: float = 1e-4
    max_kinetic_energy_fraction: float = 1e-8
    max_speed_m_per_s: float = 1e-4
    max_relative_mu_variation: float = .01


MIN_TRANSITION_CELLS = 8.0
TRANSITION_PHI_LIMIT = .9
ANGLE_ERROR_DEG = 3.0
ANGLE_FIT_SPREAD_DEG = 1.0
LAPLACE_PRESSURE_RELATIVE_TOLERANCE = .05
CHEMICAL_EQUILIBRIUM_RELATIVE_TOLERANCE = .01
