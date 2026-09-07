"""Independent equilibrium diagnostics and transient term comparison."""
import numpy as np
from .equilibrium import _integral, _measures, stationarity, chemical_stationary_form
from .p2_certification import coefficients, evaluate
from .resolution import REFERENCE_POINTS
from .free_energy import lambda_from_sigma, wall_derivative


def wall_residual(solver):
    """Strong trace L2 via exact-degree edge quadrature; sampled maximum labelled.

    Uses actual polynomial values and affine physical gradients on each boundary
    edge, independently of the weak chemical assembly. For P2, L is quartic
    along an edge, so a degree-12 rule integrates its squared norm exactly.
    """
    from .diagnostics import resolution_cell_data
    tri,v,g,*_=resolution_cell_data(solver)
    c=solver.config
    t,w=np.polynomial.legendre.leggauss((c.quadrature_degree+2)//2)
    t=(t+1)/2;w=w/2
    samples=np.linspace(0,1,65)
    normals={"bottom":[0.,-1.],"top":[0.,1.],"left":[-1.,0.],"right":[1.,0.]}
    boundaries={"bottom":(1,c.z_min),"top":(1,c.z_max),"left":(0,c.x_min),"right":(0,c.x_max)}
    co=coefficients(v) if c.phase_degree==2 else None
    integral,maximum=0.,0.
    for side in c.wetting_walls:
        axis,position=boundaries[side]
        for i,j in ((0,1),(1,2),(2,0)):
            select=np.flatnonzero(np.isclose(tri[:,i,axis],position,rtol=0,atol=1e-12)&
                                  np.isclose(tri[:,j,axis],position,rtol=0,atol=1e-12))
            if not len(select):
                continue
            length=np.linalg.norm(tri[select,j]-tri[select,i],axis=1)
            def trace(parameters):
                refs=(1-parameters[:,None])*REFERENCE_POINTS[i]+parameters[:,None]*REFERENCE_POINTS[j]
                phase=(evaluate(co[select,None,:],refs[None,:,:]) if co is not None else
                       v[select,i,None]*(1-parameters)+v[select,j,None]*parameters)
                gradient=g[select,i,None,:]*(1-parameters[None,:,None])+g[select,j,None,:]*parameters[None,:,None]
                return lambda_from_sigma(c.sigma)*c.epsilon*(gradient@normals[side])+wall_derivative(phase,c.sigma,c.theta_equilibrium_deg)
            L=trace(t)
            integral+=float(np.sum(length*np.sum(L*L*w,axis=1)))
            maximum=max(maximum,float(np.max(abs(trace(samples)))))
    from mpi4py import MPI
    return {"L2":float(np.sqrt(solver.comm.allreduce(integral,op=MPI.SUM))),
            "max_sampled":solver.comm.allreduce(maximum,op=MPI.MAX),
            "max_method":"65 points per wetting boundary edge; diagnostic, not supremum certificate",
            "quadrature_degree":c.quadrature_degree}


def state_report(solver, mu_star=None, geometry=True):
    import ufl
    from .diagnostics import Diagnostics, interface_resolution, resolution_cell_data
    row=Diagnostics(solver).measure()
    phi=solver.state.sub(2).collapse()
    mu=solver.state.sub(3).collapse()
    dx,_=_measures(solver)
    area=_integral(solver,1*dx)
    mean=_integral(solver,mu*dx)/area
    row["mu_L2_variation"]=float(np.sqrt(max(0.,_integral(solver,(mu-mean)**2*dx))))
    row["mu_coefficient_spread"]=float(np.ptp(mu.x.array))
    # An unprepared state is tested against its best L2 constant mu.
    row["stationarity"]=stationarity(solver,phi,mean if mu_star is None else mu_star)
    row["wall_residual"]=wall_residual(solver)
    row["resolution"]=interface_resolution(solver)
    ranges=resolution_cell_data(solver)[3]
    row["phi_min_certified"]=float(ranges["min"].min())
    row["phi_max_certified"]=float(ranges["max"].max())
    if geometry:
        from .benchmarks.contact_angle import observe
        row["geometry"]=observe(solver)
    return row


def transient_term_residuals(solver, mu_star, scheme="be"):
    """Separate actual transient test blocks and compare chemical sign/order.

    Block algebraic norms are debugging information only. The scientific
    stationarity gate always uses the L2 Riesz norm in state_report.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from petsc4py import PETSc
    if solver.comm.size!=1:
        raise ValueError("Term comparison currently requires one MPI rank")
    from .time_integrator import derivative_coefficients
    for coeff,value in zip(solver.coefficients,derivative_coefficients(2 if scheme=="bdf2" else 1,solver.config.dt,scheme)):
        coeff.value=value
    residual=fp.assemble_vector(fem.form(solver.F))
    residual.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
    result={}
    for index,name in enumerate(("momentum","continuity","phase","chemical")):
        _,mapping=solver.space.sub(index).collapse()
        result[name+"_algebraic_norm_diagnostic"]=float(np.linalg.norm(residual.array[mapping]))
    phi=solver.state.sub(2).collapse()
    stationary=fp.assemble_vector(fem.form(chemical_stationary_form(solver,phi,mu_star,ufl.TestFunction(phi.function_space))))
    stationary.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
    _,mapping=solver.space.sub(3).collapse()
    result["chemical_plus_stationary_algebraic_norm"]=float(np.linalg.norm(residual.array[mapping]+stationary.array))
    result["note"]="transient chemical = minus stationary EL; raw norms diagnostic only; no essential-BC lifting"
    result["quadrature_degree"]=solver.config.quadrature_degree
    result["scheme"]=scheme
    result["time_coefficients"]=[float(c.value) for c in solver.coefficients]
    stationary.destroy();residual.destroy()
    return result


class PreservationDiagnostics:
    def __init__(self,solver,phi_eq):
        import ufl
        from dolfinx import fem
        self.solver=solver
        dx,_=_measures(solver)
        phi=ufl.split(solver.state)[2]
        self.l2=fem.form((phi-phi_eq)**2*dx)
        self.h1=fem.form(((phi-phi_eq)**2+ufl.inner(ufl.grad(phi-phi_eq),ufl.grad(phi-phi_eq)))*dx)
        self.base_norm=np.sqrt(_integral(solver,phi_eq*phi_eq*dx))

    def measure(self):
        from dolfinx import fem
        from mpi4py import MPI
        l2=np.sqrt(max(0.,self.solver.comm.allreduce(float(fem.assemble_scalar(self.l2)),op=MPI.SUM)))
        h1=np.sqrt(max(0.,self.solver.comm.allreduce(float(fem.assemble_scalar(self.h1)),op=MPI.SUM)))
        return {"phi_L2_drift":float(l2),"phi_H1_drift":float(h1),
                "relative_phi_L2_drift":float(l2/max(self.base_norm,1e-30))}
