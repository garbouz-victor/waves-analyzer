"""Exact BE free-energy remainders, DIAGNOSTIC ONLY; not dissipation powers.

R(a,b) = f'(b)(b-a) - (f(b)-f(a)). The nonconvex bulk and cubic wall
remainders can be negative. Nothing is clipped or added to continuum gates.
"""
import math
from .free_energy import lambda_from_sigma


def bulk_remainder(old, new, sigma, epsilon):
    return lambda_from_sigma(sigma)/(4*epsilon)*(new-old)**2*(3*new**2+2*new*old+old**2-2)


def wall_remainder(old, new, sigma, theta_deg):
    return sigma*math.cos(math.radians(theta_deg))/4*(new-old)**2*(2*new+old)
