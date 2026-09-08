"""Immutable endpoint-subdivided BE schedules; no nonlinear input or optimizer."""
from dataclasses import dataclass, asdict
import copy
import numpy as np
from .linearized_ch import json_hash

VERSION = "parent-endpoint-phase-rate-BE-v1"


@dataclass(frozen=True)
class Interval:
    step: int
    start: float
    end: float
    dt: float
    block: int
    parent_step: int


class RateSchedule:
    def __init__(self, plan, factor=1, *, expected_hash=None, expected_count=None):
        if factor not in (1,2,4):
            raise ValueError("Only fixed 1/2/4 subdivision allowed")
        semantic=json_hash({k:v for k,v in plan.items() if k!="plan_sha256"})
        if semantic!=plan["plan_sha256"] or (expected_hash and semantic!=expected_hash):
            raise ValueError("Parent plan SHA mismatch")
        if expected_count is not None and plan["total_steps"]!=expected_count:
            raise ValueError("Parent interval count mismatch")
        parents=[]; previous=0.; ends=[]
        for block_id,block in enumerate(plan["blocks"]):
            if block["t_start"]!=previous or block["steps"]<=0 or block["dt"]<=0:
                raise ValueError("Block gap/overlap or invalid step")
            if block["t_start"]+block["steps"]*block["dt"]!=block["t_end"]:
                raise ValueError("Committed endpoint does not equal original block expression")
            for j in range(block["steps"]):
                end=block["t_start"]+(j+1)*block["dt"]
                parents.append(Interval(len(parents)+1,previous,end,block["dt"],block_id,len(parents)+1))
                previous=end
            ends.append(len(parents))
        if len(parents)!=plan["total_steps"] or previous!=plan["t_end"]:
            raise ValueError("Incorrect count/horizon")
        intervals=[]
        for parent in parents:
            a=parent.start
            for j in range(factor):
                b=parent.end if j==factor-1 else parent.start+(j+1)*(parent.end-parent.start)/factor
                intervals.append(Interval(len(intervals)+1,a,b,parent.dt/factor,parent.block,parent.step)); a=b
        for i,r in enumerate(intervals):
            if (not np.isfinite([r.start,r.end,r.dt]).all() or r.start>=r.end or
                    (i and r.start!=intervals[i-1].end)):
                raise ValueError("Invalid expanded time grid")
            if abs((r.end-r.start)-r.dt)>64*np.finfo(float).eps*max(abs(r.end),r.dt):
                raise ValueError("Clock discrepancy exceeds floating-point endpoint rounding")
        self.intervals=tuple(intervals); self.parents=tuple(parents)
        self.factor=factor; self.parent_sha256=semantic
        self.block_ends=tuple(ends); self.nsteps=len(intervals)
        self.pilot_end=min(len(parents),max(50,ends[0]+5))
        indices={0,len(parents),*ends}
        n=1
        while n<=len(parents): indices.add(n); n*=2
        indices.update(i for i in (25,26,27,28,50,self.pilot_end,*range(ends[0]+1,ends[0]+6)) if i<=len(parents))
        self.common_parent_indices=tuple(sorted(indices))
        self.field_indices=frozenset(i*factor for i in indices)
        self.checkpoint_indices=frozenset({0,self.nsteps,25*factor,50*factor,ends[0]*factor,self.pilot_end*factor} |
            set(range(200,self.nsteps+1,200)))
        self.audit_indices=frozenset({1,ends[0]*factor,ends[0]*factor+1,self.nsteps} |
            {factor*2**k for k in range(20) if factor*2**k<=self.nsteps})
        self._descriptor={"algorithm":VERSION,"parent_plan_sha256":semantic,"factor":factor,
            "level":int(np.log2(factor)),"expected_steps":self.nsteps,"final_time":previous,
            "intervals":[asdict(r) for r in intervals],"common_parent_indices":list(self.common_parent_indices),
            "checkpoint_indices":sorted(i for i in self.checkpoint_indices if i<=self.nsteps),
            "audit_indices":sorted(i for i in self.audit_indices if i<=self.nsteps),
            "clock_note":"original nominal BE dt/factor; canonical subdivided physical endpoints"}
        self.sha256=json_hash(self._descriptor)

    @property
    def descriptor(self):
        return {**copy.deepcopy(self._descriptor),"schedule_sha256":self.sha256}


def cost_projection(rows, counts, spent_s, initialization_s, limits):
    """Only measured phase-rate accepted times; never historical CHNS timing."""
    measured=np.asarray([r["timing"]["total_s"] for r in rows if r["step"]],float)
    if not len(measured) or not np.isfinite(measured).all() or np.min(measured)<=0:
        raise ValueError("Finite accepted-step timing required")
    mean,p95=float(measured.mean()),float(np.percentile(measured,95))
    per_step=1.25*max(mean,p95)
    predictions=[per_step*n+(initialization_s if n else 0.) for n in counts]
    return {"mean_s":mean,"p95_s":p95,"safety_factor":1.25,"projected_s_per_step":per_step,
        "remaining_counts":counts,"remaining_wall_s":predictions,"spent_wall_s":spent_s,
        "projected_series_wall_s":spent_s+sum(predictions),
        "series_cost_ok":spent_s+sum(predictions)<=limits["isolated_series_wall_s"]}


def scalar_convergence(integrals, changes, defects, phi_gaps, final_energies, policy):
    """Strict inherited STEP3A2 trends; no fabricated roundoff exception."""
    I,F,B=np.asarray(integrals),np.asarray(changes),np.asarray(defects)
    if len(I)!=3 or not np.isfinite(np.r_[I,F,B,phi_gaps,final_energies]).all():
        raise ValueError("Three finite complete levels required")
    gaps=abs(np.diff(I)); energy_gaps=abs(np.diff(final_energies))
    relative=[float(gaps[k]/abs(I[k+1])) if I[k+1] else None for k in (0,1)]
    closure=float(abs(B[-1]/F[-1])) if F[-1] else None
    gates={"nonzero_relaxation":bool(np.all(I>0) and np.all(F<0)),
        "dissipation_gaps_decrease":bool(gaps[1]<gaps[0]),
        "last_D_difference":relative[-1] is not None and relative[-1]<=policy["last_D_difference_relative"],
        "budget_decreases":bool(np.all(np.diff(abs(B))<0)),
        "finest_closure":closure is not None and closure<=policy["final_budget_relative"],
        "endpoint_phi_converges":bool(phi_gaps[1]<phi_gaps[0]),
        "endpoint_energy_converges":bool(energy_gaps[1]<energy_gaps[0])}
    def order(a,b): return float(np.log2(a/b)) if a>0 and b>0 else None
    return {"gates":gates,"passed":all(gates.values()),"D_absolute_gaps":gaps.tolist(),
        "D_relative_gaps":relative,"energy_gaps":energy_gaps.tolist(),"finest_budget_over_change":closure,
        "observed_orders":{"phi":order(*phi_gaps),"D_gaps":order(*gaps),
            "budget":[order(abs(B[k]),abs(B[k+1])) for k in (0,1)]}}
