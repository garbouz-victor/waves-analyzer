"""Read-only last-interval BE work decomposition; never substitutes its gate."""
import json
from pathlib import Path
import ufl
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.storage import load_checkpoint
from sloshing.multiphase.diagnostics import Diagnostics
from sloshing.multiphase.energy_validation import read_history
from sloshing.multiphase.equilibrium import energy_form,chemical_stationary_form,_integral,_measures
from sloshing.multiphase.free_energy import lambda_from_sigma
from sloshing.multiphase.be_energy_work import bulk_remainder,wall_remainder
from sloshing.multiphase.mesh import WALL_IDS
from sloshing.multiphase.provenance import run_provenance

root=Path("validation_results/step3a2/equilibrium_perturbation")
results=[]
for name in ("be_dt_1e-5","be_dt_5e-6","be_dt_2_5e-6"):
    folder=root/name
    record=json.loads((folder/"summary.json").read_text())
    c=ModelConfig(**record["config"])
    if c.time_scheme!="be" or c.rho_liquid!=c.rho_gas or c.g or c.a_x:
        raise ValueError("Derived work identity is restricted to matched-density BE, zero body force")
    s=CHNSSolver(c)
    load_checkpoint(folder/"checkpoint",s,Diagnostics(s))
    u,_,phi,_=ufl.split(s.state)
    u0,_,phi0,_=ufl.split(s.older)  # advance stores n-1 in older after acceptance
    dp,du=phi-phi0,u-u0
    dx,ds=_measures(s)
    B=(c.rho_liquid/2*ufl.inner(du,du)+lambda_from_sigma(c.sigma)*c.epsilon/2*
        ufl.inner(ufl.grad(dp),ufl.grad(dp))+bulk_remainder(phi0,phi,c.sigma,c.epsilon))*dx
    for wall in c.wetting_walls:
        B+=wall_remainder(phi0,phi,c.sigma,c.theta_equilibrium_deg)*ds(WALL_IDS[wall])
    remainder=_integral(s,B)
    work=_integral(s,c.rho_liquid*ufl.inner(du,u)*dx+chemical_stationary_form(s,phi,0.,dp))
    delta=_integral(s,c.rho_liquid/2*(ufl.inner(u,u)-ufl.inner(u0,u0))*dx+energy_form(s,phi)-energy_form(s,phi0))
    rows=read_history(folder/"history.csv")
    def power(row):
        return sum(row[key+"_dissipation"] for key in ("CH","viscous","slip"))
    current,previous=power(rows[-1]),power(rows[-2])
    weak_defect=work+c.dt*current
    gap=c.dt/2*(previous-current)
    continuum=delta+c.dt/2*(previous+current)
    results.append({"run":name,"dt":c.dt,"interval_end":s.time,"BE_remainder":remainder,
        "weak_work_defect":weak_defect,"trapezoid_minus_BE_endpoint_power":gap,
        "continuum_local_defect_from_fields":continuum,"continuum_local_defect_from_history":rows[-1]["local_budget_defect"],
        "decomposition_roundoff":continuum-(weak_defect-remainder+gap),"provenance":run_provenance(s)})
result={"scope":"LAST accepted BE interval only; read-only checkpoints, not a global/startup qualification",
    "diagnostic_only":True,"identity":"Delta E + dt*D_new + B_BE = weak_work_defect",
    "continuum_decomposition":"Delta E + trapezoid(D) = weak_work_defect - B_BE + dt/2*(D_old-D_new)",
    "B_BE":"rho/2*||du||^2 + lambda*epsilon/2*||grad dphi||^2 + integral bulk_remainder + integral wall_remainder",
    "nonnegative_remainder_not_assumed":True,"measurements":results}
with (root/"be_work_audit.json").open("x") as stream:
    json.dump(result,stream,indent=2,allow_nan=False)
print(json.dumps(result,indent=2),flush=True)
