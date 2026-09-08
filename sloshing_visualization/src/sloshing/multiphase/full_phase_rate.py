"""Validation-only full CHNS BE column reparameterization; no new physics."""
from contextlib import nullcontext
import time
import numpy as np
from .phase_rate import PhaseMetric, reconstruct_phase, snes_options

VERSION = "matched-zero-force-full-CHNS-phase-rate-BE-v1"


def require_scope(c, scheme="be"):
    if (scheme != "be" or c.rho_liquid != c.rho_gas or c.g != 0 or c.a_x != 0
            or c.phase_degree != 2 or c.quadrature_degree != 12):
        raise ValueError("Full phase-rate validation ONLY supports matched density, zero forces, P2/quad12, BE")


def physical_coefficients(old, unknown, phase_map, dt):
    result = np.array(unknown, copy=True)
    result[phase_map] = reconstruct_phase(np.asarray(old)[phase_map], result[phase_map], dt)
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite physical full candidate")
    return result


def transformed_direction(direction, phase_map, dt):
    result = np.array(direction, copy=True)
    result[phase_map] *= dt
    return result


def rate_terms(state, old, phi_old, dt, config, tags):
    """Independent transcription of the existing matched-density BE terms.

    Keys chemical_wall/navier are subterms for audits, not extra equations.
    There is no phase-state subtraction anywhere in this residual.
    """
    import ufl
    from .free_energy import bulk_derivative, lambda_from_sigma, wall_derivative
    from .material import density, density_derivative, viscosity
    from .mesh import WALL_IDS
    require_scope(config)
    c = config; mesh = state.function_space.mesh
    u, pi, rate, mu = ufl.split(state)
    u0 = ufl.split(old)[0]
    v, q, s, chi = ufl.TestFunctions(state.function_space)
    phi = phi_old+dt*rate
    dx = ufl.Measure("dx", domain=mesh, metadata={"quadrature_degree": c.quadrature_degree})
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=tags,
                     metadata={"quadrature_degree": c.quadrature_degree})
    rho, eta, drho = density(phi,c), viscosity(phi,c), density_derivative(c)
    J = -drho*c.mobility*ufl.grad(mu)  # identically zero ONLY under the scope check
    flux = rho*u+J
    momentum = (ufl.inner(rho*(u-u0)/dt,v)
        + .5*(ufl.inner(ufl.dot(ufl.grad(u),flux),v)-ufl.inner(ufl.dot(ufl.grad(v),flux),u))
        + 2*eta*ufl.inner(ufl.sym(ufl.grad(u)),ufl.sym(ufl.grad(v)))
        - pi*ufl.div(v)+phi*ufl.dot(ufl.grad(mu),v))*dx
    navier = 0*ufl.inner(u,v)*dx
    n = ufl.FacetNormal(mesh)
    for side, marker in WALL_IDS.items():
        if side not in c.no_slip_walls and side not in c.free_slip_walls:
            navier += eta/c.slip_length_m*ufl.inner(u-ufl.dot(u,n)*n,v-ufl.dot(v,n)*n)*ds(marker)
    momentum += navier
    continuity = q*ufl.div(u)*dx
    phase = (rate*s-phi*ufl.dot(u,ufl.grad(s))+c.mobility*ufl.inner(ufl.grad(mu),ufl.grad(s)))*dx
    x=ufl.SpatialCoordinate(mesh); potential=c.g*x[1]-c.a_x*x[0]
    chemical = ((mu-drho*potential-bulk_derivative(phi,c.sigma,c.epsilon))*chi
        -lambda_from_sigma(c.sigma)*c.epsilon*ufl.inner(ufl.grad(phi),ufl.grad(chi)))*dx
    wall = 0*phi*chi*dx
    for side in c.wetting_walls:
        wall -= wall_derivative(phi,c.sigma,c.theta_equilibrium_deg)*chi*ds(WALL_IDS[side])
    return {"momentum":momentum,"continuity":continuity,"phase":phase,
            "chemical":chemical+wall,"chemical_wall":wall,"navier":navier}


class CHNSPhaseRateBE:
    """Same physical layout and BCs; retained rate lives in a PRIVATE unknown."""
    def __init__(self, solver, phi0=None, *, initial_guess="semidiscrete"):
        import ufl
        from dolfinx import fem
        from dolfinx.fem.petsc import NonlinearProblem
        from petsc4py import PETSc
        require_scope(solver.config)
        if initial_guess not in ("semidiscrete","zero"):
            raise ValueError("Unknown initial rate guess")
        self.solver=solver; self.performance=None; self.space=solver.space
        if phi0 is not None and (solver.step_number!=0 or solver.time!=0.):
            raise ValueError("Explicit full initial condition only on a fresh t=0 solver")
        if phi0 is not None: solver.initialize(phi0)
        self.maps=[np.asarray(self.space.sub(i).collapse()[1]) for i in range(4)]
        self.u_map,self.pi_map,self.rate_map,self.mu_map=self.maps
        self.physical_phi_map,self.physical_mu_map=self.rate_map.copy(),self.mu_map.copy()
        self.phi_old=solver.state.sub(2).collapse()
        self.metric=PhaseMetric(solver,self.phi_old.function_space)
        self.state,self.old_state=fem.Function(self.space),fem.Function(self.space)
        self.old_state.x.array[:]=solver.state.x.array
        self.state.x.array[:]=solver.state.x.array
        mu=solver.state.x.array[self.mu_map]
        self.initial_semidiscrete_rate=self.metric.semidiscrete_rate(mu,solver.config.mobility)
        self.state.x.array[self.rate_map]=self.initial_semidiscrete_rate if initial_guess=="semidiscrete" else 0.
        self.state.x.scatter_forward(); self.old_state.x.scatter_forward()
        self.dt=fem.Constant(solver.mesh,PETSc.ScalarType(solver.config.dt))
        self.terms=rate_terms(self.state,self.old_state,self.phi_old,self.dt,solver.config,solver.tags)
        self.F=sum(self.terms[k] for k in ("momentum","continuity","phase","chemical"))
        self.J=ufl.derivative(self.F,self.state,ufl.TrialFunction(self.space))
        self.bcs=solver.bcs
        self.problem=NonlinearProblem(self.F,self.state,J=self.J,bcs=self.bcs,
            petsc_options_prefix="step3a8_full_rate_",petsc_options=snes_options(solver.config))
        self.history=[]; self.last_accepted_rate=None

    def step(self, dt, end_time=None, *, candidate_validator=None, scheme="be"):
        from dolfinx import fem
        from .nonlinear_accuracy import (synchronized_solution,check_candidate_material,
            publish_candidate,rate_arrays)
        require_scope(self.solver.config,scheme)
        if not np.isfinite(dt) or dt<=0: raise ValueError("Finite positive immutable BE dt required")
        s=self.solver; started=time.perf_counter()
        self.old_state.x.array[:]=s.state.x.array; self.old_state.x.scatter_forward()
        self.phi_old.x.array[:]=s.state.x.array[self.rate_map]; self.phi_old.x.scatter_forward()
        for mapping in (self.u_map,self.pi_map,self.mu_map):
            self.state.x.array[mapping]=s.state.x.array[mapping]
        if self.last_accepted_rate is not None:
            self.state.x.array[self.rate_map]=self.last_accepted_rate
        self.state.x.scatter_forward(); self.dt.value=dt
        with self.performance.measure("nonlinear_SNES",dt=float(dt)) if self.performance else nullcontext():
            self.problem.solve()
        snes=self.problem.solver
        if snes.getConvergedReason()<=0: raise RuntimeError("Full phase-rate SNES failed; no retry/publish")
        synchronized_solution(self.problem,self.state)
        candidate=fem.Function(s.space)
        candidate.x.array[:]=physical_coefficients(s.state.x.array,self.state.x.array,self.rate_map,dt)
        candidate.x.scatter_forward()
        target=float(end_time) if end_time is not None else s.time+dt
        if not np.isfinite(target) or target<=s.time: raise ValueError("Invalid physical endpoint")
        check_candidate_material(s,candidate)
        rate=self.state.x.array[self.rate_map].copy(); mu=self.state.x.array[self.mu_map].copy()
        arrays=rate_arrays(self.phi_old.x.array,rate,mu,dt,
            u_new=self.state.x.array[self.u_map],pi_new=self.state.x.array[self.pi_map],
            physical_state=candidate.x.array,mixed_Function=self.state.x.array,PETSc_solution=self.problem.x.array,
            rate_map=self.rate_map,mu_map=self.mu_map,u_map=self.u_map,pi_map=self.pi_map)
        result={"step":s.step_number+1,"time":target,"dt":dt,"reason":int(snes.getConvergedReason()),
            "snes_iterations":int(snes.getIterationNumber()),"residual":float(snes.getFunctionNorm()),
            "runtime_s":time.perf_counter()-started,"dt_reductions":0,"formulation":VERSION,
            "rate_L2":self.metric.norm(rate),"rate_max_abs":float(max(abs(rate))),
            "dt_times_rate_max_abs":float(dt*max(abs(rate))),"mu_L2":self.metric.norm(mu)}
        if candidate_validator is not None: candidate_validator(candidate,arrays,dt,target)
        publish_candidate(s,candidate,dt,target,"full_CHNS_phase_rate_BE")
        self.last_accepted_rate=rate; self.last_snapshot=arrays; self.history.append(result)
        return result
