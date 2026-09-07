"""Dimensional free energies, usable with scalars, NumPy or UFL expressions."""
import math


def lambda_from_sigma(sigma):
    return 3 * sigma / (2 * math.sqrt(2))


def bulk_energy(phi, sigma, epsilon):
    return lambda_from_sigma(sigma) / (4 * epsilon) * (phi**2 - 1)**2


def bulk_derivative(phi, sigma, epsilon):
    return lambda_from_sigma(sigma) / epsilon * (phi**3 - phi)


def wall_energy(phi, sigma, theta_deg):
    return -sigma * math.cos(math.radians(theta_deg)) * (3 * phi - phi**3) / 4


def wall_derivative(phi, sigma, theta_deg):
    return -3 * sigma * math.cos(math.radians(theta_deg)) * (1 - phi**2) / 4


def transition_width(epsilon):
    """Distance between phi=-0.9 and +0.9, i.e. 5–95% liquid fraction."""
    return 2 * math.sqrt(2) * math.atanh(0.9) * epsilon
