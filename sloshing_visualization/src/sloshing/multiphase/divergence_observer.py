"""New A9 acceptance policy, separate from the frozen historical A8 observer."""
import time
import numpy as np
from .full_rate_diagnostics import FullObserver, assemble_vector
from .full_rate_experiment import OneStepAudit
from .divergence_audit import DivergenceAudit
from .divergence_policy import divergence_gates


class DivergenceObserver(FullObserver):
    def __init__(self, solver, phase_metric, *, production=False):
        super().__init__(solver,phase_metric)
        self.projection = DivergenceAudit(solver)
        self.production = production

    def measure(self,row,policy,*,expanded=False):
        started = time.perf_counter()
        s = self.s; a = s.state.x.array
        velocity = a[self.maps[0]]; pressure = a[self.maps[1]]
        b = assemble_vector(self.continuity)
        vbc = float(np.max(abs(a[self.velocity_bcs]),initial=0.))
        gauge = float(np.max(abs(a[self.gauge]),initial=0.))
        result = {"F_CH":row["E_interface"]+row["E_wall"],
            "phi_L2":self.phase_metric.norm(a[self.maps[2]]),"mu_L2":self.phase_metric.norm(a[self.maps[3]]),
            "velocity_L2":self.velocity_metric.norm(velocity),"pressure_L2":self.pressure_metric.norm(pressure),
            "pressure_min_dof":float(pressure.min()),"pressure_max_dof":float(pressure.max()),
            "pressure_gauge_value":float(a[self.gauge[0]]),"velocity_BC_max":vbc,
            "velocity_BC_scaled":vbc/self.U_ref,"pressure_gauge_scaled":gauge/self.P_ref,
            "weak_continuity_algebraic":float(np.linalg.norm(b)),
            "weak_continuity_Riesz":self.pressure_metric.riesz_norm(b)}
        result.update(self.projection.measure(result["weak_continuity_Riesz"]))
        checks = {"velocity_BC":result["velocity_BC_scaled"]<=policy["BC_scaled"],
            "pressure_gauge":result["pressure_gauge_scaled"]<=policy["BC_scaled"],
            "finite_full":bool(np.isfinite(list(result.values())).all()),
            **divergence_gates(result,production=self.production)}
        result["full_physical_checks"] = {k:bool(v) for k,v in checks.items()}
        if expanded:
            adv = self.phase_metric.riesz_norm(assemble_vector(self.advection))
            diff = self.phase_metric.riesz_norm(s.config.mobility*(self.phase_metric.K@a[self.maps[3]]))
            significant = row["CH_dissipation"]>=policy["power_floor_fraction"]*policy["initial_D"]
            result["phase_coupling_audit"] = {"advection_Riesz":adv,"diffusion_Riesz":diff,
                "advective_over_diffusive":adv/diff if significant and diff>0 else None,
                "diffusion_significant":bool(significant),"diagnostic_only":True}
        result["divergence_observer_wall_s"] = time.perf_counter()-started
        return result


class DivergenceOneStepAudit(OneStepAudit):
    def __init__(self,*args,production=False,**kwargs):
        self.production = production
        super().__init__(*args,**kwargs)

    def sample(self,diag):
        if not isinstance(self.full,DivergenceObserver):
            self.full = DivergenceObserver(self.observer,self.engine.metric,production=self.production)
        return super().sample(diag)
