"""BE phase-rate change of variables, NOT a change of CH physics.

The absolute IsolatedCH and production CHNSSolver remain untouched. This
module stores rate unknowns separately and publishes only physical fields.
"""
from contextlib import nullcontext
import time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

VERSION = "isolated-CH-BE-phase-rate-v1"
ABSOLUTE_VERSION = "STEP3A4-IsolatedCH-absolute"


def reconstruct_phase(old, rate, dt):
    """Allocate a physical coefficient vector, without aliasing or corrections."""
    old, rate = np.asarray(old), np.asarray(rate)
    if old.shape != rate.shape or not np.isfinite(dt) or dt <= 0:
        raise ValueError("Compatible coefficients and finite positive BE dt required")
    result = old+dt*rate
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite reconstructed phase")
    return result


class PhaseMetric:
    """Same-quadrature sparse M/K for a guess and Riesz diagnostics; no projection."""
    def __init__(self, solver, space):
        import ufl
        from dolfinx import fem
        from dolfinx.fem import petsc as fp
        from .equilibrium import _measures
        if solver.comm.size != 1:
            raise ValueError("Validation phase-rate implementation currently requires serial layout")
        dx, _ = _measures(solver)
        trial, test = ufl.TrialFunction(space), ufl.TestFunction(space)
        matrices = []
        for form in (trial*test*dx, ufl.inner(ufl.grad(trial), ufl.grad(test))*dx):
            a = fp.assemble_matrix(fem.form(form)); a.assemble()
            ptr, col, val = a.getValuesCSR()
            matrices.append(csr_matrix((val.copy(), col.copy(), ptr.copy()), shape=a.getSize()))
            a.destroy()
        self.M, self.K = matrices
        self.lu = splu(self.M.tocsc())
        self.mass = self.M @ np.ones(self.M.shape[0])
        self.area = float(self.mass.sum())

    def norm(self, value):
        return float(np.sqrt(max(0., np.vdot(value, self.M @ value).real)))

    def riesz_norm(self, residual):
        return float(np.sqrt(max(0., np.vdot(residual, self.lu.solve(residual)).real)))

    def semidiscrete_rate(self, mu, mobility):
        return -mobility*self.lu.solve(self.K @ mu)


def snes_options(config, line_search_trace=None):
    options = {"snes_type": "newtonls", "snes_linesearch_type": "bt", "snes_stol": 0.,
        "snes_atol": config.snes_atol, "snes_rtol": config.snes_rtol, "snes_max_it": config.snes_max_it,
        "ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps",
        "snes_error_if_not_converged": True, "ksp_error_if_not_converged": True}
    if line_search_trace is not None:
        options["snes_linesearch_monitor"] = "ascii:"+str(line_search_trace)
    return options


class IsolatedCHPhaseRate:
    """Full nonlinear quartic/wetting CH in (phase_rate, mu), BE only.

First guess is the unprojected semidiscrete rate. Thereafter the previous
accepted rate is used, unscaled even at dt changes. Setting phi0=None preserves
an already initialized (including prepared constant-mu) physical observer.
"""
    def __init__(self, solver, phi0=None, *, initial_guess="semidiscrete", line_search_trace=None):
        import ufl
        from basix.ufl import element, mixed_element
        from dolfinx import fem
        from dolfinx.fem.petsc import NonlinearProblem
        from petsc4py import PETSc
        from .equilibrium import chemical_stationary_form, _measures
        c = solver.config
        if c.rho_liquid != c.rho_gas or c.g or c.a_x or c.phase_degree != 2 or c.quadrature_degree != 12:
            raise ValueError("Validation BE phase rate requires matched density, zero forcing, P2 and quadrature 12")
        if initial_guess not in ("semidiscrete", "zero"):
            raise ValueError("Unknown deterministic initial rate guess")
        self.solver, self.performance = solver, None
        self.initial_guess_policy = initial_guess
        if phi0 is not None:
            solver.initialize(phi0)  # unchanged consistent physical weak mu
            solver.step_number, solver.time = 0, 0.
            solver.current_dt, solver.current_phase = None, "initial"
        self.phi_old = solver.state.sub(2).collapse()
        self.metric = PhaseMetric(solver, self.phi_old.function_space)
        Q = element("Lagrange", solver.mesh.basix_cell(), c.phase_degree)
        self.space = fem.functionspace(solver.mesh, mixed_element([Q, Q]))
        self.state = fem.Function(self.space)
        self.dt = fem.Constant(solver.mesh, PETSc.ScalarType(c.dt))
        self.rate_map = np.asarray(self.space.sub(0).collapse()[1])
        self.mu_map = np.asarray(self.space.sub(1).collapse()[1])
        self.physical_phi_map = np.asarray(solver.space.sub(2).collapse()[1])
        self.physical_mu_map = np.asarray(solver.space.sub(3).collapse()[1])
        old_mu = solver.state.x.array[self.physical_mu_map]
        rate0 = self.metric.semidiscrete_rate(old_mu, c.mobility)
        self.state.x.array[self.rate_map] = rate0 if initial_guess == "semidiscrete" else 0.
        self.state.x.array[self.mu_map] = old_mu
        self.state.x.scatter_forward()
        self.last_accepted_rate = None
        self.initial_semidiscrete_rate = rate0.copy()
        rate, mu = ufl.split(self.state)
        test, chi = ufl.TestFunctions(self.space)
        dx, _ = _measures(solver)
        self.phi_new_expr = self.phi_old+self.dt*rate
        self.phase_form = (rate*test+c.mobility*ufl.inner(ufl.grad(mu), ufl.grad(test)))*dx
        self.chemical_form = mu*chi*dx-chemical_stationary_form(solver, self.phi_new_expr, 0., chi)
        self.F = self.phase_form+self.chemical_form
        self.J = ufl.derivative(self.F, self.state, ufl.TrialFunction(self.space))
        self.options = snes_options(c, line_search_trace)
        self.problem = NonlinearProblem(self.F, self.state, J=self.J, petsc_options_prefix="step3a5_rate_",
            petsc_options=self.options)
        self.history = []

    def rate_statistics(self, dt):
        rate = self.state.x.array[self.rate_map]
        mu = self.state.x.array[self.mu_map]
        return {"rate_L2": self.metric.norm(rate), "rate_max_abs": float(max(abs(rate))),
            "dt_times_rate_max_abs": float(dt*max(abs(rate))), "mu_L2": self.metric.norm(mu)}

    def step(self, dt, end_time=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Positive finite BE dt required")
        s = self.solver
        before = time.perf_counter()
        # phi_old is retained after this solve so its actual residual can be
        # reassembled. It is refreshed only when starting the NEXT attempt.
        self.phi_old.x.array[:] = s.state.x.array[self.physical_phi_map]
        self.phi_old.x.scatter_forward()
        if self.last_accepted_rate is not None:
            self.state.x.array[self.rate_map] = self.last_accepted_rate
        self.state.x.array[self.mu_map] = s.state.x.array[self.physical_mu_map]
        self.state.x.scatter_forward()
        self.dt.value = dt
        with self.performance.measure("nonlinear_SNES", dt=float(dt)) if self.performance else nullcontext():
            self.problem.solve()
        snes = self.problem.solver
        if snes.getConvergedReason() <= 0:
            raise RuntimeError("Phase-rate SNES failed; no publish or retry")
        self.state.x.scatter_forward()
        rate = self.state.x.array[self.rate_map].copy()
        mu = self.state.x.array[self.mu_map].copy()
        physical = reconstruct_phase(self.phi_old.x.array, rate, dt)
        if not np.isfinite(mu).all():
            raise RuntimeError("Nonfinite mu; no publish")
        # Atomic physical publication occurs only after positive convergence.
        s.older.x.array[:] = s.state.x.array
        s.state.x.array[self.physical_phi_map] = physical
        s.state.x.array[self.physical_mu_map] = mu
        s.state.x.scatter_forward(); s.older.x.scatter_forward()
        s.old.x.array[:] = s.state.x.array; s.old.x.scatter_forward()
        s._check_material()
        self.last_accepted_rate = rate
        s.step_number += 1
        s.time = float(end_time) if end_time is not None else s.time+dt
        s.current_dt, s.current_phase = dt, "isolated_ch_phase_rate_be"
        result = {"step": s.step_number, "time": s.time, "dt": dt, "reason": snes.getConvergedReason(),
            "snes_iterations": snes.getIterationNumber(), "residual": snes.getFunctionNorm(),
            "runtime_s": time.perf_counter()-before, "dt_reductions": 0,
            "formulation": VERSION, **self.rate_statistics(dt)}
        self.history.append(result)
        return result
