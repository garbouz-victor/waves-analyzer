"""STEP3A9 runner: unchanged full rate solve, independent continuity diagnostics."""
import numpy as np
from .full_rate_trajectory import FullRateTrajectory
from .divergence_observer import DivergenceObserver


class DivergenceTrajectory(FullRateTrajectory):
    def observe(self,diag):
        if self.full_observer is None:
            self.full_observer=DivergenceObserver(self.observer,self.engine.metric,production=True)
        return super().observe(diag)

    def status(self):
        result=super().status(); rows=self.journal.rows
        statistics={}
        for key in ("strong_divergence_L2","projected_divergence_L2","orthogonal_divergence_L2","projection_fraction"):
            values=np.array([r[key] for r in rows]); i=int(values.argmax())
            statistics[key]={"min":float(values.min()),"mean":float(values.mean()),
                "p50":float(np.percentile(values,50)),"p95":float(np.percentile(values,95)),
                "max":float(values[i]),"max_step":rows[i]["step"],"max_time":rows[i]["time"]}
        i=statistics["strong_divergence_L2"]["max_step"]; row=rows[i]
        result["divergence_statistics"]=statistics
        result["max_strong_divergence_state"]={k:row[k] for k in ("step","time","dt","velocity_L2",
            "E_kin","viscous_dissipation","CH_dissipation","weak_continuity_Riesz","projection_fraction")}
        result["strong_divergence_role"]="spatial diagnostic; no historical temporal envelope"
        return result
