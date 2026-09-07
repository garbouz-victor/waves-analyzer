"""Quadrature-based diagnostics, independent of any display grid.

Cumulative powers use trapezoidal TIME quadrature. For unresolved stiff startup
they can grossly overestimate dissipation; they are not certified physical heat.
The budget defect is kept visible, not corrected/clipped to force an identity.
"""
import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI
from .free_energy import bulk_energy, wall_energy, lambda_from_sigma, transition_width
from .material import density, viscosity
from .mesh import WALL_IDS


class Diagnostics:
    def __init__(self, solver):
        self.solver = solver
        self.comm = solver.comm
        c, mesh = solver.config, solver.mesh
        u, pi, phi, mu = ufl.split(solver.state)
        rho, eta = density(phi,c), viscosity(phi,c)
        x = ufl.SpatialCoordinate(mesh)
        dx = ufl.Measure("dx", domain=mesh, metadata={"quadrature_degree":c.quadrature_degree})
        ds = ufl.Measure("ds", domain=mesh, subdomain_data=solver.tags,
                         metadata={"quadrature_degree":c.quadrature_degree})
        kinetic = rho*ufl.inner(u,u)/2*dx
        gravity = rho*(c.g*x[1]-c.a_x*x[0])*dx
        interface = (bulk_energy(phi,c.sigma,c.epsilon)+lambda_from_sigma(c.sigma)*
                     c.epsilon/2*ufl.inner(ufl.grad(phi),ufl.grad(phi)))*dx
        # Multiplying zero by dx gives a domain-bearing zero form.
        wall = 0*phi*dx
        slip = 0*phi*dx
        n=ufl.FacetNormal(mesh)
        for side,marker in WALL_IDS.items():
            if side in c.wetting_walls:
                wall += wall_energy(phi,c.sigma,c.theta_equilibrium_deg)*ds(marker)
            if side not in c.no_slip_walls and side not in c.free_slip_walls:
                ut = u-ufl.dot(u,n)*n
                slip += eta/c.slip_length_m*ufl.inner(ut,ut)*ds(marker)
        forms = {"E_kin":kinetic,"E_gravity":gravity,"E_interface":interface,"E_wall":wall,
            "phase_mass":phi*dx, "liquid_volume":(1+phi)/2*dx,
            "viscous_dissipation":2*eta*ufl.inner(ufl.sym(ufl.grad(u)),ufl.sym(ufl.grad(u)))*dx,
            "CH_dissipation":c.mobility*ufl.inner(ufl.grad(mu),ufl.grad(mu))*dx,
            "slip_dissipation":slip,
            "strong_divergence_squared":ufl.div(u)**2*dx,
            "gradient_squared":ufl.inner(ufl.grad(u),ufl.grad(u))*dx,
            "wall_normal_speed_squared":ufl.dot(u,n)**2*ds,
            "omega_squared":(u[1].dx(0)-u[0].dx(1))**2*dx}
        self.forms={k:fem.form(v) for k,v in forms.items()}
        self.initial = None
        self.previous = None
        self.cumulative = {name:0.0 for name in ("viscous","CH","slip")}

    def measure(self):
        s = self.solver
        row={key:self.comm.allreduce(float(fem.assemble_scalar(form)),op=MPI.SUM)
             for key,form in self.forms.items()}
        row["time"]=s.time
        row["E_total"]=sum(row[k] for k in ("E_kin","E_gravity","E_interface","E_wall"))
        row["strong_divergence_L2"]=np.sqrt(max(0.,row.pop("strong_divergence_squared")))
        row["relative_strong_divergence"]=(row["strong_divergence_L2"]/
            max(np.sqrt(row.pop("gradient_squared")),1e-30))
        row["wall_normal_velocity_L2"]=np.sqrt(max(0.,row.pop("wall_normal_speed_squared")))
        row["omega_L2"]=np.sqrt(max(0.,row.pop("omega_squared")))
        u=s.state.sub(0).collapse()
        phi=s.state.sub(2).collapse()
        speeds=np.linalg.norm(u.x.array.reshape(-1,2),axis=1)
        row["speed_max_dof_sample"]=self.comm.allreduce(float(speeds.max()),op=MPI.MAX)
        row["phi_min_dof"]=self.comm.allreduce(float(phi.x.array.min()),op=MPI.MIN)
        row["phi_max_dof"]=self.comm.allreduce(float(phi.x.array.max()),op=MPI.MAX)
        row["wall_relaxation_dissipation"]=0.0
        if self.initial is None:
            self.initial=dict(row)
        row["mass_error_abs"]=row["phase_mass"]-self.initial["phase_mass"]
        area=(s.config.x_max-s.config.x_min)*(s.config.z_max-s.config.z_min)
        row["mass_error_relative_to_domain"]=row["mass_error_abs"]/area
        row["liquid_volume_error_relative"]=(row["liquid_volume"]-self.initial["liquid_volume"])/abs(self.initial["liquid_volume"])
        if self.previous:
            dt=row["time"]-self.previous["time"]
            for name in self.cumulative:
                self.cumulative[name]+=dt/2*(row[name+"_dissipation"]+self.previous[name+"_dissipation"])
        row.update({"cumulative_"+key+"_dissipation":v for key,v in self.cumulative.items()})
        row["energy_budget_defect"]=row["E_total"]-self.initial["E_total"]+sum(self.cumulative.values())
        row["energy_budget_relative"]=row["energy_budget_defect"]/max(abs(self.initial["E_total"]),1e-30)
        self.previous=dict(row)
        return row


def interface_resolution(solver):
    """Cell width projected along grad(phi), not diagonal h or display pixels."""
    mesh,c=solver.mesh,solver.config
    cells=np.arange(mesh.topology.index_map(2).size_local,dtype=np.int32)
    phi=ufl.split(solver.state)[2]
    center=np.array([[1/3,1/3]])
    values=fem.Expression(phi,center).eval(mesh,cells).reshape(-1)
    gradients=fem.Expression(ufl.grad(phi),center).eval(mesh,cells).reshape(-1,2)
    norms=np.linalg.norm(gradients,axis=1)
    active=(np.abs(values)<.9)&(norms>1e-12)
    tri=mesh.geometry.x[mesh.geometry.dofmap[cells],:2]
    projected=np.einsum("cij,cj->ci",tri[active],gradients[active]/norms[active,None])
    widths=np.ptp(projected,axis=1)
    hmax=solver.comm.allreduce(float(widths.max()) if widths.size else 0.,op=MPI.MAX)
    if hmax == 0:
        return {"interface_present":False,"qualified_resolution":False}
    count=transition_width(c.epsilon)/hmax
    return {"interface_present":True,"h_normal_max_m":hmax,"epsilon_over_h_normal":c.epsilon/hmax,
            "cells_across_transition_min":count,"qualified_resolution":count>=8}
