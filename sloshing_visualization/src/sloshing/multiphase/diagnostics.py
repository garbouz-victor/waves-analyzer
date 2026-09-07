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
from .energy_validation import positive_energy_scale
from .validation_policy import EnergyPolicy, MIN_TRANSITION_CELLS
from .resolution import REFERENCE_POINTS, normal_widths
from .p2_certification import polynomial_ranges, certified_normal_widths, physical_gradients, ALGORITHM_VERSION


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
            "omega_squared":(u[1].dx(0)-u[0].dx(1))**2*dx,
            "mu_integral":mu*dx,"mu_squared_integral":mu**2*dx}
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
        mu_mean=row.pop("mu_integral")/area
        row["mu_mean_Pa"]=mu_mean
        row["mu_std_Pa"]=np.sqrt(max(0.,row.pop("mu_squared_integral")/area-mu_mean**2))
        row["relative_mu_variation"]=row["mu_std_Pa"]/max(abs(mu_mean),1e-12)
        row["mass_error_relative_to_domain"]=row["mass_error_abs"]/area
        row["liquid_volume_error_relative"]=(row["liquid_volume"]-self.initial["liquid_volume"])/abs(self.initial["liquid_volume"])
        intervals={name:0. for name in self.cumulative}
        if self.previous:
            dt=row["time"]-self.previous["time"]
            for name in self.cumulative:
                intervals[name]=dt/2*(row[name+"_dissipation"]+self.previous[name+"_dissipation"])
                self.cumulative[name]+=intervals[name]
        row.update({"interval_"+key+"_dissipation":v for key,v in intervals.items()})
        row["delta_E_interval"]=row["E_total"]-self.previous["E_total"] if self.previous else 0.
        row["local_budget_defect"]=row["delta_E_interval"]+sum(intervals.values())
        row.update({"cumulative_"+key+"_dissipation":v for key,v in self.cumulative.items()})
        row["energy_budget_defect"]=row["E_total"]-self.initial["E_total"]+sum(self.cumulative.values())
        row["energy_budget_relative"]=row["energy_budget_defect"]/max(abs(self.initial["E_total"]),1e-30)
        row["E_positive_scale"]=positive_energy_scale(row)
        scale=max(positive_energy_scale(self.initial),EnergyPolicy().scale_floor)
        row["energy_budget_relative_to_initial_scale"]=row["energy_budget_defect"]/scale
        change=abs(row["E_total"]-self.initial["E_total"])
        row["energy_change_relative_applicable"]=float(change>max(EnergyPolicy().change_floor_absolute,
                                                            EnergyPolicy().change_floor_fraction*scale))
        row["energy_budget_relative_to_change"]=abs(row["energy_budget_defect"])/change if row["energy_change_relative_applicable"] else 0.
        row["energy_budget_symmetric_diagnostic"]=abs(row["energy_budget_defect"])/(change+sum(self.cumulative.values())+EnergyPolicy().change_floor_absolute)
        row["energy_budget_symmetric_diagnostic_only"]=1.0  # numeric history CSV; true
        self.previous=dict(row)
        return row


def resolution_cell_data(solver):
    """Actual reference nodal values and affine gradient samples, no raster."""
    mesh,c=solver.mesh,solver.config
    cells=np.arange(mesh.topology.index_map(2).size_local,dtype=np.int32)
    phi=solver.state.sub(2).collapse()
    nnodes=3 if c.phase_degree==1 else 6
    if not hasattr(solver,"_resolution_cell_dofs"):
        Q=phi.function_space
        native_points=Q.element.interpolation_points
        order=np.linalg.norm(REFERENCE_POINTS[:nnodes,None,:]-native_points[None,:,:],axis=-1).argmin(axis=1)
        if not np.allclose(native_points[order],REFERENCE_POINTS[:nnodes],rtol=0,atol=1e-14):
            raise ValueError("Unexpected Lagrange triangle interpolation nodes")
        solver._resolution_cell_dofs=np.array([Q.dofmap.cell_dofs(int(i))[order] for i in cells])
    values=phi.x.array[solver._resolution_cell_dofs]
    tri=mesh.geometry.x[mesh.geometry.dofmap[cells],:2]
    gradients=physical_gradients(tri,values,REFERENCE_POINTS,c.phase_degree)
    ranges=polynomial_ranges(values,c.phase_degree)
    active=(ranges["min"]<=.9)&(ranges["max"]>=-.9)
    indices=np.flatnonzero(active)
    certified=certified_normal_widths(tri[active],gradients[active,:3])
    directional=normal_widths(tri[active],gradients[active])
    return tri,values,gradients,ranges,indices,certified,directional


def interface_resolution(solver):
    """Hard gate: actual P2 ranges and certified complete-gradient-cone width."""
    c=solver.config
    tri,values,gradients,ranges,indices,certified,directional=resolution_cell_data(solver)
    widths=certified["width"]
    counts=transition_width(c.epsilon)/widths
    worst=None
    def witness(i):
        ci=int(indices[i]);triangle=tri[ci]
        result={"h_normal_m":float(widths[i]),"h_cert":float(widths[i]),
               "h_directional":float(directional[i]),"rank":solver.comm.rank,"local_cell":ci,
               "centroid":triangle.mean(axis=0).tolist(),
               "bounds":[triangle.min(axis=0).tolist(),triangle.max(axis=0).tolist()],
               "vertices":triangle.tolist(),"phi_min_exact":float(ranges["min"][ci]) if ranges["exact"][ci] else None,
               "phi_max_exact":float(ranges["max"][ci]) if ranges["exact"][ci] else None,
               "phi_min_bound":float(ranges["min"][ci]),"phi_max_bound":float(ranges["max"][ci]),
               "range_status":str(ranges["status"][ci]),
               "range_fallback_reason":str(ranges["fallback_reason"][ci]),
               "normal_width_fallback_reason":str(certified["fallback_reason"][i])}
        for side in ("min","max"):
            point=ranges[side+"_point"][ci]
            result[side+"_witness"]={"candidate_type":str(ranges[side+"_type"][ci]),
                "reference_point":point.tolist() if ranges["exact"][ci] else None,
                "physical_point":(triangle[0]+point@(triangle[1:]-triangle[0])).tolist() if ranges["exact"][ci] else None,
                "value":float(ranges[side+"_value"][ci]) if ranges["exact"][ci] else None}
        return result
    if widths.size:
        worst=witness(int(np.argmax(widths)))
    exact_indices=np.flatnonzero(ranges["exact"][indices])
    exact_worst=witness(int(exact_indices[np.argmax(widths[exact_indices])])) if len(exact_indices) else None
    all_worst=solver.comm.allgather(worst)
    all_exact_worst=solver.comm.allgather(exact_worst)
    all_counts=np.concatenate(solver.comm.allgather(counts))
    all_directional=np.concatenate(solver.comm.allgather(transition_width(c.epsilon)/directional))
    hmax=max((r["h_normal_m"] for r in all_worst if r),default=0.)
    if hmax == 0:
        return {"interface_present":False,"qualified_resolution":False,"active_cell_count":0,
                "active_cell_count_exact":0,"algorithm_version":ALGORITHM_VERSION}
    count=transition_width(c.epsilon)/hmax
    worst=max((r for r in all_worst if r),key=lambda r:r["h_normal_m"])
    quantiles=np.percentile(all_counts,[5,50,95])
    result={"interface_present":True,"h_normal_max_m":hmax,"epsilon_over_h_normal":c.epsilon/hmax,
            "cells_across_transition_min":count,"qualified_resolution":count>=MIN_TRANSITION_CELLS,
            "cells_across_transition_p05":float(quantiles[0]),"cells_across_transition_p50":float(quantiles[1]),
            "cells_across_transition_p95":float(quantiles[2]),"active_cell_count":len(all_counts),
            "active_cell_count_exact":solver.comm.allreduce(len(exact_indices),op=MPI.SUM),
            "active_cell_count_conservative_fallback":solver.comm.allreduce(int((~ranges["exact"][indices]).sum()),op=MPI.SUM),
            "diameter_fallback_active_cells":solver.comm.allreduce(int(certified["diameter_fallback"].sum()),op=MPI.SUM),
            "diameter_fallback_local_cells_by_rank":solver.comm.allgather(indices[certified["diameter_fallback"]].tolist()),
            "worst_cell":worst,"worst_certified_cell":worst,
            "worst_exact_active_cell":max((r for r in all_exact_worst if r),key=lambda r:r["h_cert"],default=None),
            "time":solver.time,"algorithm_version":ALGORITHM_VERSION,
            "activation":"CG1 exact vertices; P2 vertices + edge extrema + interior stationary point; labelled Bernstein fallback",
            "h_normal_method":"maximum over complete affine-gradient cone; diameter for origin/near-origin hull",
            "hard_gate_metric":"cells_across_transition_certified_min",
            "cells_across_transition_certified_min":count,
            "cells_across_transition_directional_min":float(all_directional.min())}
    for name,data in (("certified",all_counts),("directional",all_directional)):
        for pct,value in zip(("p05","p50","p95"),np.percentile(data,[5,50,95])):
            result[f"cells_across_transition_{name}_{pct}"]=float(value)
    return result
