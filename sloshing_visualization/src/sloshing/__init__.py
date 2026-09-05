"""Numerical physics first; renderers consume saved solutions."""

from .config import SimulationConfig
from .solver import SloshingSolver

__all__ = ["SimulationConfig", "SloshingSolver"]
