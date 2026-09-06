"""First-order Lagrangian displacement at FIXED reference points."""

import numpy as np


def linear_lagrangian_displacement(times, velocity_at_seeds):
    """Trapezoid is exact for the chosen piecewise-linear temporal interpolant.

    For samples of a smooth analytic function its integration error is O(dt²).
    Input shape (time, particle, xy); no Eulerian v2 is invented.
    """
    times = np.asarray(times)
    v = np.asarray(velocity_at_seeds)
    if v.ndim != 3 or v.shape[-1] != 2 or len(v) != len(times) or np.any(np.diff(times) <= 0):
        raise ValueError("Expected increasing times and (time,particle,2) velocities")
    increments = np.diff(times)[:, None, None]*(v[1:]+v[:-1])/2
    return np.concatenate((np.zeros_like(v[:1]), np.cumsum(increments, axis=0)))


def choose_magnification(displacement, target=.075):
    q = float(np.percentile(np.linalg.norm(displacement, axis=-1), 95))
    return (target/q if q > 0 else 1.), q
