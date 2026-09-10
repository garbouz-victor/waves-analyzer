"""All dimensional parameters are SI; energies are per unit density and span."""

from dataclasses import asdict, dataclass
import hashlib
import json
import math


MESH_PRESETS = {"coarse": (40, 70), "medium": (80, 140), "fine": (120, 210)}
VISCOSITY_PRESETS = (0.001, 0.01, 0.1, 1.0)


@dataclass(frozen=True)
class SimulationConfig:
    a: float = 1.0
    d: float = 10.0
    g: float = 9.81
    alpha_deg: float = 2.0
    nu: float = 0.01
    t_end: float = 5.0
    dt: float = 0.0025
    integrator: str = "midpoint"
    snapshot_dt: float = 0.025
    mesh: str = "medium"
    nx: int = None
    nz: int = None
    x_grading: float = 1.8
    z_grading: float = 4.0
    visualization_nx: int = 81
    visualization_nz: int = 241
    small_slope_limit: float = 0.3
    weak_divergence_tolerance: float = 1e-9
    wall_tolerance: float = 1e-12
    volume_tolerance: float = 1e-10
    contact_tolerance: float = 1e-12
    energy_relative_tolerance: float = 1e-8
    divergence_relative_warning: float = 0.15

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, (int, float)) and not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        for name in ("a", "d", "g", "nu", "dt", "snapshot_dt", "small_slope_limit",
                     "weak_divergence_tolerance", "wall_tolerance", "volume_tolerance",
                     "contact_tolerance", "energy_relative_tolerance", "divergence_relative_warning"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.t_end < 0 or abs(self.alpha_deg) >= 90:
            raise ValueError("t_end must be nonnegative and |alpha_deg| < 90 degrees")
        if self.mesh not in MESH_PRESETS:
            raise ValueError(f"mesh must be one of {tuple(MESH_PRESETS)}")
        if self.integrator not in ("midpoint", "sdirk2"):
            raise ValueError("integrator must be midpoint or sdirk2")
        for name in ("nx", "nz", "visualization_nx", "visualization_nz"):
            n = getattr(self, name)
            if n is not None and (not isinstance(n, int) or isinstance(n, bool) or n < 4):
                raise ValueError(f"{name} must be an integer >= 4")
        if not (0 <= self.x_grading <= 5 and 0 <= self.z_grading <= 8):
            raise ValueError("grading outside supported range: x in [0,5], z in [0,8]")
        for name in ("t_end", "snapshot_dt"):
            steps = getattr(self, name) / self.dt
            if not math.isclose(steps, round(steps), rel_tol=0, abs_tol=1e-8):
                raise ValueError(f"{name} must be an integer multiple of dt")
        if self.snapshot_dt < self.dt:
            raise ValueError("snapshot_dt must be >= dt")

    @property
    def resolution(self):
        nx, nz = MESH_PRESETS[self.mesh]
        return self.nx or nx, self.nz or nz

    @property
    def slope(self):
        return math.tan(math.radians(self.alpha_deg))

    @property
    def nsteps(self):
        return round(self.t_end / self.dt)

    @property
    def save_every(self):
        return round(self.snapshot_dt / self.dt)

    def to_json(self):
        return json.dumps(asdict(self), sort_keys=True)

    @classmethod
    def from_json(cls, path):
        with open(path, encoding="utf-8") as handle:
            return cls(**json.load(handle))

    @property
    def run_name(self):
        digest = hashlib.sha256(self.to_json().encode()).hexdigest()[:10]
        nx, nz = self.resolution
        return f"run_nu{self.nu:g}_alpha{self.alpha_deg:g}_{nx}x{nz}_{digest}"
