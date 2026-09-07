"""Explicit dimensional parameters for STEP 3A; no water/air default fiction."""
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class ModelConfig:
    name: str
    rho_liquid: float
    rho_gas: float
    mu_liquid: float
    mu_gas: float
    sigma: float
    epsilon: float
    mobility: float
    slip_length_m: float
    theta_equilibrium_deg: float
    g: float
    a_x: float
    length_scale_m: float
    velocity_scale_m_per_s: float
    x_min: float
    x_max: float
    z_min: float
    z_max: float
    nx: int
    nz: int
    phase_degree: int
    dt: float
    t_end: float
    x_grading: float = 0.0
    z_grading: float = 0.0
    wetting_walls: tuple = ("left", "right", "bottom")
    no_slip_walls: tuple = ("bottom",)
    free_slip_walls: tuple = ("top",)
    time_scheme: str = "bdf2"
    quadrature_degree: int = 12
    snes_atol: float = 1e-10
    snes_rtol: float = 1e-9
    snes_max_it: int = 30
    max_unknowns: int = 100000
    refinement_levels: int = 0
    refinement_circle_x: float = 0.0
    refinement_circle_z: float = 0.0
    refinement_circle_radius: float = 0.25
    refinement_band_m: float = 0.12

    def __post_init__(self):
        for key,value in asdict(self).items():
            if isinstance(value,(int,float)) and not math.isfinite(value):
                raise ValueError(f"{key} must be finite")
        for key in ("nx","nz","phase_degree","quadrature_degree","refinement_levels","max_unknowns"):
            if not isinstance(getattr(self,key),int):
                raise ValueError(f"{key} must be an integer")
        positive = ("rho_liquid", "rho_gas", "mu_liquid", "mu_gas", "sigma",
                    "epsilon", "mobility", "slip_length_m", "length_scale_m",
                    "velocity_scale_m_per_s", "dt", "t_end")
        for name in positive:
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.rho_liquid < self.rho_gas or self.g < 0:
            raise ValueError("This convention requires rho_liquid >= rho_gas and g >= 0")
        if not 0 < self.theta_equilibrium_deg < 180:
            raise ValueError("Young angle must lie strictly between 0 and 180 degrees")
        if self.x_min >= self.x_max or self.z_min >= self.z_max:
            raise ValueError("Invalid rectangle")
        if min(self.nx, self.nz) < 2 or self.phase_degree not in (1, 2):
            raise ValueError("Need >=2 intervals and phase degree 1 or 2")
        if min(self.x_grading, self.z_grading) < 0:
            raise ValueError("Grading must be nonnegative")
        if self.refinement_levels not in (0,1,2,3) or min(self.refinement_circle_radius,self.refinement_band_m)<=0:
            raise ValueError("Invalid local interface-refinement request")
        if self.time_scheme not in ("be", "bdf2"):
            raise ValueError("Supported schemes: be, bdf2 (BE startup)")
        if self.quadrature_degree < 6 * self.phase_degree:
            # Quartic energy requires 4p; chemical residual/Jacobian with tests
            # and material/capillary products benefit from the stricter 6p rule.
            raise ValueError("Use quadrature_degree >= 6*phase_degree")
        sides = {"left", "right", "bottom", "top"}
        for name in ("wetting_walls", "no_slip_walls", "free_slip_walls"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
            if not set(getattr(self, name)) <= sides:
                raise ValueError(f"Unknown boundary in {name}")
        if set(self.no_slip_walls) & set(self.free_slip_walls):
            raise ValueError("A wall cannot be both no-slip and free-slip")

    @classmethod
    def from_json(cls, path):
        return cls(**json.loads(Path(path).read_text()))

    def changed(self, **kwargs):
        return replace(self, **kwargs)

    def as_dict(self):
        return asdict(self)

    def fingerprint(self):
        return hashlib.sha256(json.dumps(self.as_dict(), sort_keys=True).encode()).hexdigest()

    def groups(self):
        L, U, sigma = self.length_scale_m, self.velocity_scale_m_per_s, self.sigma
        dr = self.rho_liquid - self.rho_gas
        lc = math.sqrt(sigma / (dr * self.g)) if dr * self.g > 0 else None
        return {"L_m": L, "U_m_per_s": U,
                "Re": self.rho_liquid * U * L / self.mu_liquid,
                "Fr": U / math.sqrt(self.g * L) if self.g else None,
                "We": self.rho_liquid * U**2 * L / sigma,
                "Bo": dr * self.g * L**2 / sigma,
                "Ca": self.mu_liquid * U / sigma, "Cn": self.epsilon / L,
                "slip_ratio": self.slip_length_m / L,
                "density_ratio": self.rho_liquid / self.rho_gas,
                "viscosity_ratio": self.mu_liquid / self.mu_gas,
                "Pe_CH": U * L**2 / (self.mobility * sigma),
                "Pe_CH_convention": "mu_scale=sigma/L; U L^2/(M sigma)",
                "capillary_length_m": lc,
                "capillary_length_infinite": lc is None,
                "half_width_over_lc": (self.x_max-self.x_min)/2/lc if lc else 0.0,
                "epsilon_over_lc": self.epsilon/lc if lc else 0.0}
