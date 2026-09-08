"""Read-only full residual/Jacobian bridge audit on arbitrary tiny FE fields."""
import numpy as np
from .full_phase_rate import CHNSPhaseRateBE,physical_coefficients,transformed_direction
from .full_rate_diagnostics import assemble_vector,constrained_maps
from .full_rate_policy import POLICY


def sparse_matrix(form,bcs=()):
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from scipy.sparse import csr_matrix
    a=fp.assemble_matrix(fem.form(form),bcs=list(bcs)); a.assemble()
    ptr,col,val=a.getValuesCSR()
    result=csr_matrix((val.copy(),col.copy(),ptr.copy()),shape=a.getSize()); a.destroy()
    return result


def relative(a,b): return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))


def audit(solver):
    import ufl
    from dolfinx import fem
    from .weak_form import residual
    from .free_energy import wall_derivative
    from .material import viscosity
    from .mesh import WALL_IDS
    e=CHNSPhaseRateBE(solver,lambda x:.2+.05*np.cos(2*x[0])*np.cos(3*x[1]),initial_guess="zero")
    rng=np.random.default_rng(38217); dt=.03; e.dt.value=dt
    e.old_state.x.array[e.u_map]=.01*rng.normal(size=len(e.u_map))
    e.state.x.array[e.u_map]=.02*rng.normal(size=len(e.u_map))
    e.state.x.array[e.pi_map]=.03*rng.normal(size=len(e.pi_map))
    e.state.x.array[e.rate_map]=.1*rng.normal(size=len(e.rate_map))
    e.state.x.array[e.mu_map]=.1*rng.normal(size=len(e.mu_map))
    e.state.x.scatter_forward(); e.old_state.x.scatter_forward()
    physical=fem.Function(e.space)
    physical.x.array[:]=physical_coefficients(e.old_state.x.array,e.state.x.array,e.rate_map,dt)
    physical.x.scatter_forward()
    oldF=residual(physical,e.old_state,e.old_state,(1/dt,-1/dt,0.),solver.config,solver.tags)
    oldR=assemble_vector(oldF); rateR=assemble_vector(e.F)
    errors={k:relative(rateR[m],oldR[m]) for k,m in zip(("momentum","continuity","phase","chemical"),e.maps)}
    u,pi,phi,mu=ufl.split(physical); v,q,s,chi=ufl.TestFunctions(e.space)
    ds=ufl.Measure("ds",domain=solver.mesh,subdomain_data=solver.tags,metadata={"quadrature_degree":12})
    wall=sum(-wall_derivative(phi,solver.config.sigma,solver.config.theta_equilibrium_deg)*chi*ds(WALL_IDS[w])
             for w in solver.config.wetting_walls)
    n=ufl.FacetNormal(solver.mesh); eta=viscosity(phi,solver.config)
    navier=sum(eta/solver.config.slip_length_m*ufl.inner(u-ufl.dot(u,n)*n,v-ufl.dot(v,n)*n)*ds(marker)
        for w,marker in WALL_IDS.items() if w not in solver.config.no_slip_walls and w not in solver.config.free_slip_walls)
    for name,form in (("chemical_wall",wall),("navier",navier)):
        errors[name]=relative(assemble_vector(e.terms[name]),assemble_vector(form))
    oldJ=sparse_matrix(ufl.derivative(oldF,physical,ufl.TrialFunction(e.space)))
    newJ=sparse_matrix(e.J)
    dy=rng.normal(size=len(e.state.x.array)); dx=transformed_direction(dy,e.rate_map,dt)
    action=newJ@dy; oldaction=oldJ@dx
    jac_error=relative(action,oldaction)
    base=e.state.x.array.copy(); fd=[]
    for h in (1e-2,5e-3,2.5e-3,1.25e-3,1e-5):
        e.state.x.array[:]=base+h*dy; e.state.x.scatter_forward(); plus=assemble_vector(e.F)
        e.state.x.array[:]=base-h*dy; e.state.x.scatter_forward(); minus=assemble_vector(e.F)
        fd.append({"h":h,"relative_error":relative((plus-minus)/(2*h),action)})
    e.state.x.array[:]=base; e.state.x.scatter_forward()
    bc_u,bc_p=constrained_maps(solver)
    constrained_rate=np.intersect1d(np.r_[bc_u,bc_p],e.rate_map)
    fd_trend=all(b["relative_error"]<=a["relative_error"]*.4+1e-12 for a,b in zip(fd[:3],fd[1:4]))
    gates={"term_residuals":max(errors.values())<=POLICY["bridge_relative"],
        "Jacobian_column_transformation":jac_error<=POLICY["bridge_relative"],
        "FD_second_order_until_roundoff":fd_trend,
        "FD_final":fd[-1]["relative_error"]<=POLICY["FD_final_relative"],
        "BC_identical_objects":e.bcs is solver.bcs,
        "one_pressure_gauge":len(bc_p)==1,"no_rate_BC":len(constrained_rate)==0}
    return {"terms_relative":errors,"Jacobian_action_relative":jac_error,"finite_difference":fd,
        "BC_velocity_count":len(bc_u),"pressure_gauge_dof":int(bc_p[0]),
        "gates":gates,"passed":all(gates.values()),"arbitrary_nonzero_fields":True}
