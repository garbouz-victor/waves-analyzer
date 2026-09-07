"""STEP3A4 fixed scientific limits and role-specific execution budgets."""
from dataclasses import asdict, dataclass
import time
from .linearized_ch import json_hash


LINEAR_LIMITS = {"endpoint_relative_L2_error": 1e-3,
    "max_checkpoint_relative_L2_error": 1e-3,
    "max_quadratic_energy_relative_error": 1e-3,
    "integrated_D2_relative_error": 1e-2}


@dataclass(frozen=True)
class ValidationPolicy:
    analysis_stage_wall_s: float = 3600.
    analysis_total_wall_s: float = 7200.
    optimizer_cpu_s: float = 1800.
    optimizer_evaluations: int = 5000
    optimizer_rounds: int = 24
    optimizer_max_splits: int = 3
    isolated_level0_wall_s: float = 7200.
    isolated_series_wall_s: float = 28800.
    isolated_pilot_steps: int = 50
    coupling_wall_s: float = 7200.
    linear_mass_domain: float = 1e-11
    nonlinear_mass_domain: float = 1e-10
    certified_cells: float = 8.
    energy_growth_J_per_m: float = 1e-9
    weak_work_J_per_m: float = 1e-12
    final_budget_relative: float = .05
    last_D_difference_relative: float = .05
    coupling_phi_relative: float = 1e-3
    coupling_energy_relative: float = 1e-3
    coupling_power_relative: float = 1e-2
    coupling_kinetic_over_excess: float = 1e-4
    coupling_hydro_over_CH: float = 1e-4
    coupling_power_floor_fraction: float = 1e-8
    coupling_field_floor_fraction: float = 1e-12
    t_end: float = 1e-4

    def record(self):
        return {**asdict(self), "linear_limits": dict(LINEAR_LIMITS),
            "version": "step3a4-policy-v1", "nonlinear_SNESTolerances_changed": False}

    @property
    def fingerprint(self):
        return json_hash(self.record())


def linear_gates(errors, full_sparse=False, policy=ValidationPolicy()):
    gates = {key: 0 <= errors.get(key, float("inf")) <= limit
             for key, limit in LINEAR_LIMITS.items()}
    if full_sparse:
        gates["mass_domain"] = 0 <= errors.get("max_mass_domain", float("inf")) <= policy.linear_mass_domain
    return gates


def authorizations(*, reduced=False, sparse=False, isolated_cost_ok=False,
                   isolated_series_pass=False, probe_cost_ok=False, coupling_pass=False):
    """Scientific roles are distinct; there is deliberately no generic nonlinear flag."""
    return {"reduced_linear_verified": bool(reduced),
        "full_sparse_linear_verified": bool(reduced and sparse),
        "isolated_CH_execution_authorized": bool(reduced and sparse and isolated_cost_ok),
        "full_CHNS_probe_authorized": bool(reduced and sparse and isolated_series_pass and probe_cost_ok),
        "full_model_temporal_qualification": "not_qualified",  # limited conditional bridge is not this
        "conditional_CH_dominated_transfer": bool(reduced and sparse and isolated_series_pass and coupling_pass)}


class AnalysisBudget:
    def __init__(self, previously_used_s=0., policy=ValidationPolicy()):
        self.started = time.perf_counter()
        self.used = previously_used_s
        self.policy = policy

    def check(self, projected_remaining_s=0.):
        elapsed = time.perf_counter()-self.started
        if (elapsed+projected_remaining_s > self.policy.analysis_stage_wall_s or
                self.used+elapsed+projected_remaining_s > self.policy.analysis_total_wall_s):
            raise RuntimeError("STOP: independent linear analysis cost gate")


def run_sparse_analysis(verifier, *, reduced_pass, analysis_budget):
    """Regression boundary: full-CHNS cost is not an input to linear execution."""
    if not reduced_pass:
        raise ValueError("Reduced linear verification required before full sparse")
    analysis_budget.check()
    result = verifier()
    analysis_budget.check()
    return result
