"""Mechanical pressure jump checked against independently fitted interface radius."""
import json
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI
from ..initialization import droplet
from ..material import density_derivative
from ..interface import contour_segments, fit_circle, connected_polylines
from .common import run_case
from ..energy_validation import read_history, validate_energy
from ..validation_policy import LAPLACE_PRESSURE_RELATIVE_TOLERANCE, CHEMICAL_EQUILIBRIUM_RELATIVE_TOLERANCE


def measure(solver):
    c=solver.config
    segments=contour_segments(solver)
    components=connected_polylines(segments)
    if len(components)!=1:
        raise RuntimeError(f"Laplace benchmark requires one droplet; found {len(components)} components")
    circle=fit_circle(segments)
    u,pi,phi,mu=ufl.split(solver.state)
    x=ufl.SpatialCoordinate(solver.mesh)
    psi=c.g*x[1]-c.a_x*x[0]
    mu0=mu-density_derivative(c)*psi
    p=pi+phi*mu0-(c.rho_liquid+c.rho_gas)/2*psi
    dx=ufl.Measure("dx",domain=solver.mesh,metadata={"quadrature_degree":c.quadrature_degree})
    def integral(expr):
        return solver.comm.allreduce(float(fem.assemble_scalar(fem.form(expr*dx))),op=MPI.SUM)
    radius=ufl.sqrt((x[0]-circle["center_x"])**2+(x[1]-circle["center_z"])**2)
    inside=ufl.conditional(ufl.lt(radius,circle["radius"]-4*c.epsilon),1.,0.)
    outside=ufl.conditional(ufl.gt(radius,circle["radius"]+4*c.epsilon),1.,0.)
    vin,vout=integral(inside),integral(outside)
    if min(vin,vout)<=0:
        raise RuntimeError("No pure-phase averaging region >=4 epsilon from interface")
    pin,pout=integral(p*inside)/vin,integral(p*outside)/vout
    area=(c.x_max-c.x_min)*(c.z_max-c.z_min)
    mu_average=integral(mu)/area
    chemical_std=float(np.sqrt(max(0.,integral((mu-mu_average)**2)/area)))
    expected=c.sigma/circle["radius"]
    return {**circle,"pressure_inside_Pa":pin,"pressure_outside_Pa":pout,
            "pressure_jump_Pa":pin-pout,"sigma_over_R_Pa":expected,
            "pressure_relative_error":abs(pin-pout-expected)/expected,
            "mu_std_Pa":chemical_std,"mu_mean_Pa":mu_average,
            "relative_mu_std":chemical_std/max(abs(mu_average),1e-30)},components


def run(config,output,radius=None):
    r=radius or config.refinement_circle_radius
    c=config.changed(refinement_circle_radius=r)
    solver,rows,result=run_case(c,droplet(c.epsilon,r),output)
    return solver,summarize(solver,result,output)


def summarize(solver,result,output):
    measurement,components=measure(solver)
    energy_gate=validate_energy(read_history(Path(output)/"history.csv"))
    checks={"pressure_error_below_5_percent":measurement["pressure_relative_error"]<LAPLACE_PRESSURE_RELATIVE_TOLERANCE,
            "near_chemical_equilibrium":measurement["relative_mu_std"]<CHEMICAL_EQUILIBRIUM_RELATIVE_TOLERANCE,
            "interface_resolved":result.get("interface_resolved_all_times",False),
            "mass_conservation":energy_gate["mass_ok"],
            "no_unexplained_energy_growth":energy_gate["energy_growth_ok"],
            "energy_budget_closure":energy_gate["energy_budget_ok"]}
    status="passed" if all(checks.values()) else "failed"
    summary={"status":status,"qualification_status":status,"checks":checks,"measurement":measurement,
             "energy_validation":energy_gate,
             "scope":"one radius/epsilon, not a convergence study","run":result}
    if solver.comm.rank==0:
        folder=Path(output)
        (folder/"laplace.json").write_text(json.dumps(summary,indent=2)+"\n")
        (folder/"interface.json").write_text(json.dumps([p.tolist() for p in components])+"\n")
        print(json.dumps({"laplace":measurement,"checks":checks},indent=2),flush=True)
    return summary


def analyze_complete(output):
    from ..config import ModelConfig
    from ..solver import CHNSSolver
    from ..diagnostics import Diagnostics
    from ..storage import load_checkpoint
    historical=Path(__file__).resolve().parents[4]/"validation_results/step3"
    if Path(output).resolve()==historical.resolve() or historical.resolve() in Path(output).resolve().parents:
        raise ValueError("Historical STEP 3 files are read-only; write a separate step3a1 audit")
    result=json.loads((Path(output)/"summary.json").read_text())
    if result["status"]!="complete":
        raise ValueError("Only a complete solver checkpoint can be analyzed")
    solver=CHNSSolver(ModelConfig(**result["config"]))
    load_checkpoint(Path(output)/"checkpoint",solver,Diagnostics(solver))
    return summarize(solver,result,output)
