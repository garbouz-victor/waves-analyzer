"""Measure prepared-state BE/BDF2 algebraic rounding floor before tolerance policy."""
import json
from pathlib import Path
import numpy as np
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.initialization import flat_interface
from sloshing.multiphase.equilibrium import solve_equilibrium
from sloshing.multiphase.equilibrium_diagnostics import transient_term_residuals
from sloshing.multiphase.provenance import run_provenance

c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(nx=8,nz=8,dt=1e-5,t_end=4e-5)
s=CHNSSolver(c)
s.initialize(flat_interface(c.epsilon))
p=solve_equilibrium(s)
s.initialize_prepared_equilibrium(p)
rows=[transient_term_residuals(s,p.metadata["mu_star"],scheme) for scheme in ("be","bdf2")]
result={"scope":"tiny prepared flat state; investigate over-tight 1e-13 transient SNES atol",
        "provenance":run_provenance(s),"stationarity":p.metadata["stationarity"],"blocks":rows,
        "roundoff_scale_per_unit_mass_weight":float(8*np.finfo(float).eps*4/c.dt),
        "analysis":"BDF2 forms a weighted sum of O(phi/dt) terms; floating arithmetic need not cancel identical fields. Existing transient default atol=1e-10 exceeds this tiny-mesh algebraic rounding floor. Scientific field/mass/energy thresholds are unchanged; stationary augmented atol remains 1e-13."}
output=Path("validation_results/step3a2/conditioning")
output.mkdir(parents=True,exist_ok=True)
with (output/"tiny_bdf2_roundoff.json").open("x") as stream:
    json.dump(result,stream,indent=2)
print(json.dumps(result,indent=2),flush=True)
