"""Full-CHNS validation specialization of the accepted STEP3A7 journal runner."""
import copy
import time
import numpy as np
from .phase_rate_trajectory import RateTrajectory,ScientificFailure,CostStop,physical_gates
from .phase_rate_history import archive,atomic_json,read_json,file_hash,validate_complete
from .linearized_ch import array_hash
from .full_rate_diagnostics import FullObserver


class FullRateTrajectory(RateTrajectory):
    def __init__(self,*args,full_policy,**kwargs):
        self.full_policy=full_policy; self.full_observer=None
        # The isolated expanded audit has the wrong block layout and omits advection.
        super().__init__(*args,expanded_audit=False,**kwargs)

    def initialize(self):
        if self.s.step_number!=0 or self.s.time!=0.: raise ValueError("Fresh full level starts at t=0")
        self.update_observer(self.s.state.x.array,self.s.state.x.array,0.,0,None)
        row,res,timing=self.observe(self.diag)
        gates=physical_gates(row,res,None,self.policy,self.target_mass,self.area,self.expected_crossings)
        if not all(gates.values()): raise ScientificFailure("Initial full physical gates: "+str(gates))
        self.journal.append({**row,"step":0,"dt":None,"block":None,"solver":None,"work":None,
            "physical_gates":gates,"resolution":res,"timing":timing,
            "benchmark":"FULL nonlinear matched-density zero-force CHNS phase-rate BE"})
        self.save_field(0); self.checkpoint()

    def observe(self,diag):
        row,resolution,timing=super().observe(diag)
        started=time.perf_counter()
        if self.full_observer is None:
            self.full_observer=FullObserver(self.observer,self.engine.metric)
        row.update(self.full_observer.measure(row,self.full_policy,
            expanded=self.observer.step_number in self.schedule.audit_indices))
        timing["full_constraint_audit_s"]=time.perf_counter()-started
        if not all(row["full_physical_checks"].values()):
            raise ScientificFailure("Full candidate constraints failed: "+str(row))
        return row,resolution,timing

    def physical_arrays(self):
        a=super().physical_arrays(); e=self.engine
        a.update(u=self.s.state.x.array[e.u_map].copy(),pi=self.s.state.x.array[e.pi_map].copy(),
            u_map=e.u_map.copy(),pi_map=e.pi_map.copy())
        return a

    def save_field(self,step):
        a=self.physical_arrays(); dt=self.s.current_dt
        old_phi=self.s.older.x.array[self.engine.physical_phi_map].copy()
        physical={k:a[k] for k in ("u","pi","phi","mu","phase_rate")}
        physical["phi_old"]=old_phi
        if step:
            physical.update(increment_from_rate=dt*a["phase_rate"],
                increment_from_physical_states=a["phi"]-old_phi,
                phase_rate_recovered_from_rounded_phi=(a["phi"]-old_phi)/dt)
        return archive(self.folder/"fields"/f"{step:06d}",physical,
            {"identity":self.identity,"step":step,"time":self.s.time,"dt":dt,
             "rate_role":"retained full BE unknown" if step else "semidiscrete initial guess; no accepted rate yet"})

    def restore(self,arrays,meta):
        super().restore(arrays,meta)
        for key in ("u_map","pi_map"):
            if not np.array_equal(arrays[key],getattr(self.engine,key)):
                raise ValueError("Full checkpoint DOF layout mismatch")
            m=getattr(self.engine,key); self.engine.state.x.array[m]=self.s.state.x.array[m]
        self.engine.state.x.scatter_forward()

    def run_until(self,count,*,wall_remaining_s=float("inf")):
        if not self.s.step_number<=count<=self.schedule.nsteps: raise ValueError("Invalid full prefix")
        started=time.perf_counter()
        while self.s.step_number<count:
            if time.perf_counter()-started>=wall_remaining_s:
                self.checkpoint(); raise CostStop("Full-CHNS frozen wall budget exhausted")
            row=self.advance()
            if row["step"]%10==0 or row["step"]==count:
                print(f"full step={row['step']}/{self.schedule.nsteps} t={row['time']:.10g} "
                    f"Newton={row['solver']['snes_iterations']} D={row['CH_dissipation']:.9g} "
                    f"Ekin={row['E_kin']:.3g} mass={row['mass_error_relative_to_domain']:.3g} "
                    f"cells={row['cells_across_transition_certified_min']:.9g}",flush=True)
        self.checkpoint()
        if count==self.schedule.nsteps:
            self.source_guard()
            if not all(all(r["full_physical_checks"].values()) for r in self.journal.rows):
                raise ScientificFailure("Full physical gate missing at completion")
            marker={"identity":self.identity,"accepted_steps":count,"expected_steps":count,
                "final_time":self.s.time,"final_state_sha256":array_hash(self.s.state.x.array),
                "scalar_history_sha256":file_hash(self.journal.path)}
            validate_complete(self.folder,self.identity,self.schedule,marker=marker)
            atomic_json(self.folder/"COMPLETE.json",marker)
        return self.status()

    def status(self):
        result=super().status(); rows=self.journal.rows; last=rows[-1]
        result.update(full_physical_status="passed" if all(all(r["full_physical_checks"].values()) for r in rows) else "failed",
            final_F_CH=last["F_CH"],integrated_D_visc=last["cumulative_viscous_dissipation"],
            integrated_D_slip=last["cumulative_slip_dissipation"],
            max_weak_continuity=max(r["weak_continuity_Riesz"] for r in rows),
            max_strong_divergence=max(r["strong_divergence_L2"] for r in rows),
            max_relative_divergence=max(r["relative_strong_divergence"] for r in rows),
            max_speed=max(r["speed_max_dof_sample"] for r in rows))
        return result
