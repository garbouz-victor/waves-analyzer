"""Analytic benchmark initial fields, not production pre-equilibration."""
import numpy as np


def flat_interface(epsilon, height=0.0, alpha_deg=0.0):
    angle = np.deg2rad(alpha_deg)
    # Positive alpha: liquid interface is HIGHER on the RIGHT.
    return lambda x: np.tanh((height+x[0]*np.sin(angle)-x[1]*np.cos(angle)) /
                             (np.sqrt(2)*epsilon))


def droplet(epsilon, radius, center=(0.0, 0.0)):
    return lambda x: np.tanh((radius-np.sqrt((x[0]-center[0])**2+(x[1]-center[1])**2)) /
                             (np.sqrt(2)*epsilon))


def sessile_drop(epsilon, radius, theta_deg, wall_z=0.0):
    return droplet(epsilon, radius, (0.0, wall_z-radius*np.cos(np.deg2rad(theta_deg))))


def effective_acceleration(alpha_deg, g):
    return g*np.tan(np.deg2rad(alpha_deg))
