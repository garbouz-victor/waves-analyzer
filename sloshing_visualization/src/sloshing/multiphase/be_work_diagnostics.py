"""Signed BE work remainder versus physical continuum power, diagnostic only."""


class BEWorkDiagnostics:
    def __init__(self, solver):
        import ufl
        from dolfinx import fem
        from .equilibrium import _measures, energy_form, chemical_stationary_form
        from .free_energy import lambda_from_sigma
        from .be_energy_work import bulk_remainder, wall_remainder
        from .mesh import WALL_IDS
        self.solver = solver
        c = solver.config
        if c.rho_liquid != c.rho_gas or c.g or c.a_x:
            raise ValueError("BE work identity restricted to matched density, zero body force")
        u, _, phi, _ = ufl.split(solver.state)
        old_u, _, old_phi, _ = ufl.split(solver.older)
        du, dp = u-old_u, phi-old_phi
        dx, ds = _measures(solver)
        B = (c.rho_liquid/2*ufl.inner(du, du)+lambda_from_sigma(c.sigma)*c.epsilon/2*
             ufl.inner(ufl.grad(dp), ufl.grad(dp))+bulk_remainder(old_phi, phi, c.sigma, c.epsilon))*dx
        for wall in c.wetting_walls:
            B += wall_remainder(old_phi, phi, c.sigma, c.theta_equilibrium_deg)*ds(WALL_IDS[wall])
        work = c.rho_liquid*ufl.inner(du, u)*dx+chemical_stationary_form(solver, phi, 0., dp)
        change = c.rho_liquid/2*(ufl.inner(u, u)-ufl.inner(old_u, old_u))*dx
        change += energy_form(solver, phi)-energy_form(solver, old_phi)
        self.forms = {"B_BE": fem.form(B), "energy_work": fem.form(work), "Delta_E": fem.form(change)}

    def measure(self, dt, old, new):
        from dolfinx import fem
        from mpi4py import MPI
        result = {k: self.solver.comm.allreduce(float(fem.assemble_scalar(v)), op=MPI.SUM) for k, v in self.forms.items()}
        def power(row):
            return sum(row[name+"_dissipation"] for name in ("CH", "viscous", "slip"))
        old_power, new_power = power(old), power(new)
        trap = dt/2*(old_power+new_power)
        endpoint = dt*new_power
        weak = result["energy_work"]+endpoint
        continuum = result["Delta_E"]+trap
        gap = trap-endpoint
        result.update(time=self.solver.time, step=self.solver.step_number, dt=dt,
            D_old=old_power, D_new=new_power, trapezoid_integral=trap, BE_endpoint_integral=endpoint,
            weak_work_defect=weak, continuum_local_defect=continuum,
            trapezoid_gap=gap, decomposition_roundoff=continuum-(weak-result["B_BE"]+gap),
            weak_work_identity_roundoff=result["Delta_E"]+endpoint+result["B_BE"]-weak,
            diagnostic_only=True, remainder_label="signed discrete BE work remainder; NOT physical heat")
        return result
