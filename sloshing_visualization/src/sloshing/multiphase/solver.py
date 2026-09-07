"""Benchmark CHNS engine; DOLFINx/SNES, never the old scikit-fem solver."""
import time
import numpy as np
import ufl
import basix
from basix.ufl import element, mixed_element
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem, LinearProblem
from mpi4py import MPI
from petsc4py import PETSc
from .mesh import make_mesh, WALL_IDS
from .boundary import velocity_pressure_bcs
from .weak_form import residual
from .time_integrator import derivative_coefficients, IntegrationSchedule
from .free_energy import bulk_derivative, wall_derivative, lambda_from_sigma
from .material import density_derivative, assert_admissible


class CHNSSolver:
    def __init__(self, config, comm=MPI.COMM_WORLD):
        self.config, self.comm = config, comm
        self.schedule=IntegrationSchedule(config)
        self.mesh, self.tags, self.facets, self.mesh_report = make_mesh(config, comm)
        cell = self.mesh.basix_cell()
        self.space = fem.functionspace(self.mesh, mixed_element([
            element("Lagrange", cell, 2, shape=(2,)), element("Lagrange", cell, 1),
            element("Lagrange", cell, config.phase_degree),
            element("Lagrange", cell, config.phase_degree)]))
        self.state, self.old, self.older = [fem.Function(self.space) for _ in range(3)]
        self.coefficients = [fem.Constant(self.mesh, PETSc.ScalarType(0)) for _ in range(3)]
        self.bcs = velocity_pressure_bcs(self.space, config, self.facets)
        self.F = residual(self.state, self.old, self.older, self.coefficients, config, self.tags)
        self.problem = NonlinearProblem(self.F, self.state, bcs=self.bcs,
            petsc_options_prefix="step3_", petsc_options={
                "snes_type": "newtonls", "snes_linesearch_type": "bt",
                "snes_atol": config.snes_atol, "snes_rtol": config.snes_rtol,
                "snes_stol": 0.0, "snes_max_it": config.snes_max_it,
                "ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps",
                "snes_error_if_not_converged": True, "ksp_error_if_not_converged": True})
        self.step_number, self.time = 0, 0.0
        self.current_phase="initial"
        self.current_dt=None
        self.history_spacing=None
        self.logs = []
        points, _ = basix.make_quadrature(self.mesh.basix_cell(), config.quadrature_degree)
        self.phase_quadrature = fem.Expression(ufl.split(self.state)[2], points)

    def initialize(self, phase, velocity=None):
        self.state.x.array[:] = 0
        self.state.sub(2).interpolate(phase)
        if velocity is not None:
            self.state.sub(0).interpolate(velocity)
        self.state.x.scatter_forward()
        # Consistent weak chemical potential, including natural wall variation.
        Q, mapping = self.space.sub(3).collapse()
        trial, test = ufl.TrialFunction(Q), ufl.TestFunction(Q)
        phi = ufl.split(self.state)[2]
        cfg = self.config
        x = ufl.SpatialCoordinate(self.mesh)
        dx = ufl.Measure("dx", domain=self.mesh, metadata={"quadrature_degree": cfg.quadrature_degree})
        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.tags,
                         metadata={"quadrature_degree": cfg.quadrature_degree})
        rhs = ((bulk_derivative(phi, cfg.sigma, cfg.epsilon)
                + density_derivative(cfg)*(cfg.g*x[1]-cfg.a_x*x[0]))*test
               + lambda_from_sigma(cfg.sigma)*cfg.epsilon*ufl.dot(ufl.grad(phi), ufl.grad(test))) * dx
        for wall in cfg.wetting_walls:
            rhs += wall_derivative(phi, cfg.sigma, cfg.theta_equilibrium_deg)*test*ds(WALL_IDS[wall])
        projection = LinearProblem(trial*test*dx, rhs, petsc_options_prefix="step3_mu_init_",
            petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"})
        mu = projection.solve()
        self.state.x.array[mapping] = mu.x.array
        self.state.x.scatter_forward()
        for previous in (self.old, self.older):
            previous.x.array[:] = self.state.x.array
            previous.x.scatter_forward()
        self._check_material()

    def _check_material(self):
        ncells = self.mesh.topology.index_map(2).size_local
        values = self.phase_quadrature.eval(self.mesh, np.arange(ncells, dtype=np.int32))
        failure = 0
        try:
            assert_admissible(values, self.config)
        except RuntimeError:
            failure = 1
        if self.comm.allreduce(failure, op=MPI.MAX):
            raise RuntimeError("Inadmissible material coefficient / nonfinite phi; no clipping applied")

    def initialize_prepared_equilibrium(self, prepared):
        """Preserve certified constant mu exactly; no chemical re-projection."""
        from .equilibrium import validate_prepared
        validate_prepared(self,prepared)
        _,mapping=self.space.sub(2).collapse()
        if len(mapping)!=len(prepared.phi_coefficients):
            raise ValueError("Prepared equilibrium coefficient layout mismatch")
        self.state.x.array[:]=0.
        self.state.x.array[mapping]=prepared.phi_coefficients
        _,mu_mapping=self.space.sub(3).collapse()
        self.state.x.array[mu_mapping]=prepared.metadata["mu_star"]
        self.state.x.scatter_forward()
        for previous in (self.old,self.older):
            previous.x.array[:]=self.state.x.array
            previous.x.scatter_forward()
        self.step_number,self.time=0,0.
        self.current_phase,self.current_dt,self.history_spacing="initial",None,None
        self.equilibrium_source_fingerprint=prepared.metadata["fingerprint"]
        self._check_material()

    def advance(self):
        started = time.perf_counter()
        number = self.step_number+1
        interval=self.schedule.step(self.step_number)
        if interval["scheme"]=="bdf2" and (self.history_spacing is None or
                not np.isclose(self.history_spacing,interval["dt"],rtol=1e-12,atol=1e-15)):
            raise RuntimeError("Refusing constant-step BDF2 with unequal history spacing")
        for c, value in zip(self.coefficients, derivative_coefficients(
                2 if interval["scheme"]=="bdf2" else 1,interval["dt"],interval["scheme"])):
            c.value = value
        self.problem.solve()
        snes = self.problem.solver
        if snes.getConvergedReason() <= 0:
            raise RuntimeError(f"SNES failed: {snes.getConvergedReason()}")
        self.state.x.scatter_forward()
        self._check_material()
        self.step_number, self.time = number, interval["time"]
        self.current_phase,self.current_dt=interval["phase"],interval["dt"]
        self.history_spacing=interval["dt"]
        row = {"step": number, "time": self.time,
               "newton_iterations": snes.getIterationNumber(),
               "linear_iterations": snes.getLinearSolveIterations(),
               "residual": snes.getFunctionNorm(), "reason": snes.getConvergedReason(),
               "runtime_s": time.perf_counter()-started, "dt_reductions": 0,
               "dt":self.current_dt,"scheme_phase":self.current_phase}
        self.logs.append(row)
        self.older.x.array[:] = self.old.x.array
        self.old.x.array[:] = self.state.x.array
        self.older.x.scatter_forward()
        self.old.x.scatter_forward()
        return row
