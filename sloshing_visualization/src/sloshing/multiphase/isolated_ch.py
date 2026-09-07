"""Isolated CH temporal BENCHMARK, u=0; never a production CHNS replacement."""
import time
from contextlib import nullcontext
import numpy as np


class IsolatedCH:
    def __init__(self, solver, phi0):
        import ufl
        from basix.ufl import element, mixed_element
        from dolfinx import fem
        from dolfinx.fem.petsc import NonlinearProblem
        from petsc4py import PETSc
        from .equilibrium import chemical_stationary_form, _measures
        self.solver = solver
        self.performance = None
        c = solver.config
        if c.g or c.a_x or c.rho_liquid != c.rho_gas:
            raise ValueError("Isolated CH benchmark requires matched density, zero body force")
        Q = element("Lagrange", solver.mesh.basix_cell(), c.phase_degree)
        space = fem.functionspace(solver.mesh, mixed_element([Q, Q]))
        self.state, self.old = fem.Function(space), fem.Function(space)
        self.dt = fem.Constant(solver.mesh, PETSc.ScalarType(c.dt))
        phi, mu = ufl.split(self.state)
        old_phi, _ = ufl.split(self.old)
        test, chi = ufl.TestFunctions(space)
        dx, _ = _measures(solver)
        self.phase_form = ((phi-old_phi)/self.dt*test+c.mobility*ufl.inner(ufl.grad(mu), ufl.grad(test)))*dx
        self.chemical_form = mu*chi*dx-chemical_stationary_form(solver, phi, 0., chi)
        self.F = self.phase_form+self.chemical_form
        self.problem = NonlinearProblem(self.F, self.state, petsc_options_prefix="step3a3_isolated_",
            petsc_options={"snes_type": "newtonls", "snes_linesearch_type": "bt", "snes_stol": 0.,
                "snes_atol": c.snes_atol, "snes_rtol": c.snes_rtol, "snes_max_it": c.snes_max_it,
                "ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps",
                "snes_error_if_not_converged": True, "ksp_error_if_not_converged": True})
        # Existing production weak mu initialization, same space/quadrature/wall.
        solver.initialize(phi0)
        self.state.sub(0).interpolate(solver.state.sub(2).collapse())
        self.state.sub(1).interpolate(solver.state.sub(3).collapse())
        self.state.x.scatter_forward()
        self.old.x.array[:] = self.state.x.array
        self.old.x.scatter_forward()
        solver.step_number, solver.time = 0, 0.

    def step(self, dt, end_time=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Positive dt required")
        s = self.solver
        before = time.perf_counter()
        self.dt.value = dt
        with self.performance.measure("nonlinear_SNES", dt=float(dt)) if self.performance else nullcontext():
            self.problem.solve()
        snes = self.problem.solver
        if snes.getConvergedReason() <= 0:
            raise RuntimeError("Isolated CH SNES failed; no timestep retry")
        self.state.x.scatter_forward()
        # Publish ONLY into the benchmark observer. Production advance is never
        # called; this solver is explicitly labelled isolated CH in every run.
        s.older.x.array[:] = s.state.x.array
        s.state.sub(2).interpolate(self.state.sub(0).collapse())
        s.state.sub(3).interpolate(self.state.sub(1).collapse())
        s.state.x.scatter_forward(); s.older.x.scatter_forward()
        s.old.x.array[:] = s.state.x.array; s.old.x.scatter_forward()
        s._check_material()
        self.old.x.array[:] = self.state.x.array; self.old.x.scatter_forward()
        s.time = float(end_time) if end_time is not None else s.time+dt
        s.step_number += 1
        s.current_dt, s.current_phase = dt, "isolated_ch_be"
        return {"step": s.step_number, "time": s.time, "dt": dt,
                "snes_iterations": snes.getIterationNumber(), "residual": snes.getFunctionNorm(),
                "runtime_s": time.perf_counter()-before, "dt_reductions": 0,
                "benchmark": "isolated CH temporal benchmark; u=0"}
