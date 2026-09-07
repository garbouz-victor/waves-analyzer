"""Affine density is required for exact AGG density/phase mass compatibility."""
import numpy as np


def density(phi, config):
    return (config.rho_liquid + config.rho_gas) / 2 + density_derivative(config) * phi


def density_derivative(config):
    return (config.rho_liquid - config.rho_gas) / 2


def viscosity(phi, config):
    return (config.mu_liquid + config.mu_gas) / 2 + (config.mu_liquid-config.mu_gas)/2 * phi


def assert_admissible(phi_values, config):
    if not np.all(np.isfinite(phi_values)):
        raise RuntimeError("Non-finite phase field")
    if np.min(density(phi_values, config)) <= 0 or np.min(viscosity(phi_values, config)) <= 0:
        raise RuntimeError("Nonpositive material coefficient from phase overshoot; no clipping applied")
