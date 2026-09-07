"""No timestep: independent semidiscrete CH/free-energy identity at perturbed t=0."""
import json
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.equilibrium import load_prepared,zero_mass_perturbations,energy_form,_integral,_measures
from sloshing.multiphase.provenance import run_provenance

folder=Path("validation_results/step3a2/equilibrium60/prepared_symmetric")
meta=json.loads((folder/"prepared.json").read_text())
s=CHNSSolver(ModelConfig(**meta["config"]))
p=load_prepared(folder,s,meta["target_mass"])
s.initialize_prepared_equilibrium(p)
phi_eq=s.state.sub(2).collapse()
psi,_=zero_mass_perturbations(s,phi_eq)["smooth_bulk"]
perturbed=fem.Function(phi_eq.function_space)
perturbed.x.array[:]=phi_eq.x.array+1e-3*psi.x.array
perturbed.x.scatter_forward()
s.initialize(perturbed)
phi,mu=s.state.sub(2).collapse(),s.state.sub(3).collapse()
Q=phi.function_space
dx,_=_measures(s)
trial,test=ufl.TrialFunction(Q),ufl.TestFunction(Q)
rate_problem=LinearProblem(trial*test*dx,-s.config.mobility*ufl.dot(ufl.grad(mu),ufl.grad(test))*dx,
    petsc_options_prefix="step3a2_initial_rate_",petsc_options={"ksp_type":"preonly","pc_type":"lu",
    "pc_factor_mat_solver_type":"mumps","ksp_error_if_not_converged":True})
rate=rate_problem.solve()
power=_integral(s,s.config.mobility*ufl.inner(ufl.grad(mu),ufl.grad(mu))*dx)
mu_rate=_integral(s,mu*rate*dx)
energy_rate=_integral(s,ufl.derivative(energy_form(s,phi),phi,rate))
excess=_integral(s,energy_form(s,phi))-_integral(s,energy_form(s,phi_eq))
result={"scope":"perturbed t=0 semidiscrete diagnostic; no transient advance, not a replacement energy gate",
    "equilibrium_source_fingerprint":meta["fingerprint"],"provenance":run_provenance(s),
    "D_CH":power,"mu_times_phase_rate":mu_rate,"free_energy_directional_rate":energy_rate,
    "phase_rate_L2":float(np.sqrt(_integral(s,rate*rate*dx))),"phase_mass_rate":_integral(s,rate*dx),
    "chemical_identity_defect":energy_rate-mu_rate,"CH_identity_defect":mu_rate+power,
    "semidiscrete_energy_identity_defect":energy_rate+power,
    "relative_semidiscrete_defect":abs(energy_rate+power)/power,
    "initial_excess_energy":excess,"excess_energy_over_initial_power_s":excess/power,
    "definition":"(phi_dot_h,s)=-M(grad mu,grad s); independently differentiate F_h in direction phi_dot_h",
    "diagnostic_only":True}
output=Path("validation_results/step3a2/equilibrium_perturbation/initial_rate_audit.json")
with output.open("x") as stream:
    json.dump(result,stream,indent=2,allow_nan=False)
print(json.dumps(result,indent=2),flush=True)
