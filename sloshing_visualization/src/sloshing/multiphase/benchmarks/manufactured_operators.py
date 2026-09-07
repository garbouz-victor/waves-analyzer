"""Independent analytic coupled-operator consistency, not a solution-error MMS.

SymPy supplies time-dependent unmatched-density fields and sources. Measure the
L2 Riesz representatives of FEM weak defects; raw algebraic residual norms would
vanish with mesh size even for wrong PDE signs and are NOT used as evidence.
"""
import json
from pathlib import Path
import sympy as sp
import numpy as np
import ufl
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from ..solver import CHNSSolver
from ..weak_form import residual


def analytic(config):
    x,z,t=sp.symbols("x z t",real=True)
    decay=sp.exp(-t)
    psi=sp.Rational(1,100)*decay*sp.sin(sp.pi*x)**2*sp.sin(sp.pi*z)**2
    vel=sp.Matrix([sp.diff(psi,z),-sp.diff(psi,x)])
    phi=sp.Rational(3,20)*decay*sp.cos(sp.pi*x)*sp.cos(sp.pi*z)
    pressure=sp.Rational(1,5)*decay*sp.sin(sp.pi*x)*sp.sin(sp.pi*z)
    lam=3*config.sigma/(2*sp.sqrt(2))
    mu=lam/config.epsilon*(phi**3-phi)-lam*config.epsilon*(sp.diff(phi,x,2)+sp.diff(phi,z,2))
    rho=(config.rho_liquid+config.rho_gas)/2+(config.rho_liquid-config.rho_gas)/2*phi
    rho_prime=(config.rho_liquid-config.rho_gas)/2
    eta=(config.mu_liquid+config.mu_gas)/2+(config.mu_liquid-config.mu_gas)/2*phi
    grad=lambda value:sp.Matrix([sp.diff(value,x),sp.diff(value,z)])
    jac=vel.jacobian([x,z])
    stress=eta*(jac+jac.T)
    J=-rho_prime*config.mobility*grad(mu)
    source_phi=sp.diff(phi,t)+vel.dot(grad(phi))-config.mobility*(sp.diff(mu,x,2)+sp.diff(mu,z,2))
    divstress=sp.Matrix([sp.diff(stress[i,0],x)+sp.diff(stress[i,1],z) for i in range(2)])
    # The manufactured CH source introduces a density source rho_prime*Sphi.
    # The skew strong momentum therefore includes half this source times u.
    source_v=(rho*sp.diff(vel,t)+jac*(rho*vel+J)+rho_prime*source_phi*vel/2
              -divstress+grad(pressure)+phi*grad(mu))
    assert sp.simplify(sp.diff(vel[0],x)+sp.diff(vel[1],z))==0
    fields=[vel[0],vel[1],pressure,phi,mu]
    numerical=[sp.lambdify((x,z,t),v,"numpy",cse=True) for v in fields]
    mapping={"sin":ufl.sin,"cos":ufl.cos,"exp":ufl.exp,"sqrt":ufl.sqrt}
    loads=[sp.lambdify((x,z,t),v,modules=[mapping],cse=True) for v in [*source_v,source_phi]]
    return numerical,loads


def measure(config,n,wrong_pressure_sign=False):
    c=config.changed(nx=n,nz=n,x_min=0,x_max=1,z_min=0,z_max=1,
        refinement_levels=0,rho_liquid=1.2,rho_gas=1.,mu_liquid=.2,mu_gas=.1,
        epsilon=.1,sigma=.02,mobility=.01,wetting_walls=(),free_slip_walls=(),
        no_slip_walls=("left","right","bottom","top"),g=0.,a_x=0.,dt=1e-5)
    solver=CHNSSolver(c)
    fields,loads=analytic(c)
    time=.2
    for function,t in zip((solver.state,solver.old,solver.older),(time,time-c.dt,time-2*c.dt)):
        function.sub(0).interpolate(lambda p:np.vstack([f(p[0],p[1],t) for f in fields[:2]]))
        for i,func in enumerate(fields[2:],start=1):
            function.sub(i).interpolate(lambda p,f=func:f(p[0],p[1],t))
        function.x.scatter_forward()
    for constant,value in zip(solver.coefficients,(1.5/c.dt,-2/c.dt,.5/c.dt)):
        constant.value=value
    F=residual(solver.state,solver.old,solver.older,solver.coefficients,c,solver.tags)
    v,q,s,chi=ufl.TestFunctions(solver.space)
    x=ufl.SpatialCoordinate(solver.mesh)
    force=[f(x[0],x[1],time) for f in loads]
    dx=ufl.Measure("dx",domain=solver.mesh,metadata={"quadrature_degree":16})
    F-=ufl.inner(ufl.as_vector(force[:2]),v)*dx+force[2]*s*dx
    if wrong_pressure_sign:
        F+=2*ufl.split(solver.state)[1]*ufl.div(v)*dx
    dv,dp,df,dm=ufl.TrialFunctions(solver.space)
    mass=(ufl.inner(dv,v)+dp*q+df*s+dm*chi)*dx
    problem=LinearProblem(mass,F,bcs=solver.bcs,petsc_options_prefix="step3_mms_riesz_",
        petsc_options={"ksp_type":"preonly","pc_type":"lu","pc_factor_mat_solver_type":"mumps"})
    representative=problem.solve()
    result={"n":n,"h":1/n,"wrong_pressure_sign":wrong_pressure_sign}
    for name,field in zip(("momentum","divergence","CH","chemical"),ufl.split(representative)):
        square=solver.comm.allreduce(float(fem.assemble_scalar(fem.form(ufl.inner(field,field)*dx))),op=MPI.SUM)
        result[name+"_defect_L2"]=float(np.sqrt(max(0.,square)))
    # This symmetric manufactured velocity has a near-exact weak P1 divergence
    # constraint after interpolation. Do not demand monotonic decay of roundoff.
    strong=fem.assemble_scalar(fem.form(ufl.div(ufl.split(solver.state)[0])**2*dx))
    result["strong_divergence_L2"]=float(np.sqrt(solver.comm.allreduce(float(strong),op=MPI.SUM)))
    return result


def study(config,output):
    rows=[measure(config,n) for n in (6,12,24)]
    wrong=measure(config,24,wrong_pressure_sign=True)
    checks={"all_operator_defects_decrease":all(
        b[key]<a[key] for a,b in zip(rows,rows[1:])
        for key in ("momentum_defect_L2","strong_divergence_L2","CH_defect_L2","chemical_defect_L2")),
        "weak_divergence_roundoff":max(r["divergence_defect_L2"] for r in rows)<1e-12,
        "pressure_sign_negative_control_detected":wrong["momentum_defect_L2"]>3*rows[-1]["momentum_defect_L2"]}
    result={"status":"passed" if all(checks.values()) else "failed","checks":checks,
            "rows":rows,"negative_control":wrong,
            "scope":"coupled spatial operator consistency; NOT a full solution-error convergence MMS",
            "density_ratio":1.2,"viscosity_ratio":2.,"forcing_in_normal_runs":False}
    if MPI.COMM_WORLD.rank==0:
        folder=Path(output);folder.mkdir(parents=True,exist_ok=True)
        (folder/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps(result,indent=2),flush=True)
    if result["status"]!="passed":
        raise RuntimeError("Manufactured-operator consistency failed")
    return result
