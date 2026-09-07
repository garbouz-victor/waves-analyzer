"""Deterministic analysis-only BE schedule planning on a reduced linear model.

No nonlinear observed state is accepted by this API. Small dense matrices are
Krylov projections (or synthetic tests), never the full production CH operator.
The final production plan additionally requires full sparse BE verification.
"""
from dataclasses import asdict, dataclass
import copy
import numpy as np
from scipy.linalg import eig, expm, lu_factor, lu_solve
from .linearized_ch import array_hash, json_hash

ALGORITHM_VERSION = "fixed-block-linear-BE-planner-v1"


@dataclass(frozen=True)
class PlannerPolicy:
    z_target: float = .05
    vanished_power_fraction: float = 1e-4
    endpoint_relative_tolerance: float = 1e-3
    checkpoint_relative_tolerance: float = 1e-3
    energy_relative_tolerance: float = 1e-3
    integrated_power_relative_tolerance: float = 1e-2
    max_candidate_refinements: int = 12
    max_analysis_steps: int = 50000


class ReducedCHReference:
    def __init__(self, rate, energy_matrix, power_matrix, seed, mass_matrix=None):
        self.A, self.H, self.D = (np.asarray(a, dtype=float) for a in (rate, energy_matrix, power_matrix))
        self.seed = np.asarray(seed, dtype=float)
        self.size = len(self.seed)
        self.M = np.eye(self.size) if mass_matrix is None else np.asarray(mass_matrix)
        if self.size > 4096 or any(a.shape != (self.size, self.size) for a in (self.A, self.H, self.D, self.M)):
            raise ValueError("Only compatible SMALL reduced matrices are allowed")
        if not all(np.isfinite(a).all() for a in (self.A, self.H, self.D, self.M, self.seed)):
            raise ValueError("Nonfinite reduced model")
        self.values, self.vectors = eig(self.A)
        self.weights = np.linalg.solve(self.vectors, self.seed)
        self.condition = float(np.linalg.cond(self.vectors))
        if self.condition > 1e10:
            raise ValueError("Ill-conditioned modal exponential; require direct expm reference")
        if max(self.values.real) > 1e-8*max(1., max(abs(self.values))):
            raise ValueError("Unstable reduced spectrum")
        self.fingerprint = json_hash({k: array_hash(v) for k, v in
            {"rate": self.A, "energy": self.H, "power": self.D, "mass": self.M, "seed": self.seed}.items()})
        self._lu = {}

    def norm(self, q):
        return float(np.sqrt(max(0., np.vdot(q, self.M @ q).real)))

    def energy(self, q):
        return float(.5*np.vdot(q, self.H @ q).real)

    def power(self, q):
        return float(np.vdot(q, self.D @ q).real)

    def state(self, t):
        q = self.vectors @ (np.exp(t*self.values)*self.weights)
        if np.linalg.norm(q.imag) > 1e-10*max(np.linalg.norm(q.real), 1e-300):
            raise ValueError("Complex reference state")
        return q.real

    def verify_modal_exponential(self, times):
        return max(self.norm(self.state(t)-expm(t*self.A) @ self.seed)/
                   max(self.norm(self.state(t)), 1e-14*self.norm(self.seed)) for t in times)

    def integrated_power(self, horizon):
        """Analytical integral of D2(t), NOT inferred from the energy change."""
        matrix = self.vectors.conj().T @ self.D @ self.vectors
        products = self.weights.conj()[:, None]*matrix*self.weights[None, :]
        rates = self.values.conj()[:, None]+self.values[None, :]
        factors = np.full(rates.shape, complex(horizon))
        np.divide(np.expm1(horizon*rates), rates, out=factors, where=abs(rates) > 1e-14)
        return float(np.sum(products*factors).real)

    def be_step(self, q, dt):
        if dt not in self._lu:
            self._lu[dt] = lu_factor(np.eye(self.size)-dt*self.A)
        return lu_solve(self._lu[dt], q)


def refine_blocks(blocks, factor):
    if factor not in (1, 2, 4) and (factor < 1 or factor & (factor-1)):
        raise ValueError("Plan refinement must be a power of two")
    return [{**b, "dt": b["dt"]/factor, "steps": b["steps"]*factor} for b in blocks]


def candidate_blocks(reference, horizon, rho_relevant, policy):
    if horizon <= 0 or rho_relevant <= 0:
        raise ValueError("Positive horizon and relevant radius required")
    dt0 = policy.z_target/rho_relevant
    amplitudes = abs(reference.weights)**2*np.real(np.diag(reference.vectors.conj().T @ reference.D @ reference.vectors))
    power0 = reference.power(reference.seed)
    t, previous = 0., dt0
    blocks = []
    for step in range(policy.max_analysis_steps):
        if t >= horizon:
            break
        active = amplitudes*np.exp(2*t*reference.values.real) >= policy.vanished_power_fraction*power0
        rho = max(abs(reference.values.real[active]), default=0.)
        wanted = policy.z_target/rho if rho else horizon/100
        wanted = min(wanted, horizon/100)
        dt = dt0 if step == 0 else min(2*previous, max(dt0, wanted))
        # Powers of two produce reusable sparse factors and exact refinement.
        dt = dt0*2**max(0, int(np.floor(np.log2(max(dt/dt0, 1.)))))
        dt = min(dt, horizon-t)
        if dt <= 0 or t+dt == t:
            raise RuntimeError("Unrepresentable planned interval")
        end = horizon if horizon-(t+dt) <= 1e-14*horizon else t+dt
        if blocks and dt == blocks[-1]["dt"]:
            blocks[-1]["steps"] += 1; blocks[-1]["t_end"] = end
        else:
            blocks.append({"block": len(blocks), "t_start": t, "t_end": end,
                           "dt": dt, "steps": 1, "active_stiffness": float(rho)})
        t, previous = end, dt
    else:
        raise RuntimeError("Linear proposal exceeded finite analysis-step guard")
    return blocks


def verify_schedule(reference, blocks):
    q = reference.seed.copy()
    t = 0.; integral = 0.; old_power = reference.power(q)
    errors, energy_errors, block_errors = [], [], []
    for block in blocks:
        local = []
        for step in range(block["steps"]):
            dt = block["dt"]
            t = block["t_start"]+(step+1)*dt
            q = reference.be_step(q, dt)
            exact = reference.state(t)
            error = reference.norm(q-exact)/max(reference.norm(exact), 1e-14*reference.norm(reference.seed))
            energy_error = abs(reference.energy(q)-reference.energy(exact))/max(abs(reference.energy(exact)), 1e-14*abs(reference.energy(reference.seed)))
            power = reference.power(q)
            integral += dt*(old_power+power)/2
            old_power = power
            errors.append(error); energy_errors.append(energy_error); local.append(error)
        block_errors.append(float(max(local)))
    exact_integral = reference.integrated_power(blocks[-1]["t_end"])
    reference._lu.clear()
    return {"endpoint_relative_L2_error": float(errors[-1]),
            "max_checkpoint_relative_L2_error": float(max(errors)),
            "max_quadratic_energy_relative_error": float(max(energy_errors)),
            "integrated_D2_relative_error": float(abs(integral/exact_integral-1)),
            "linear_BE_integrated_D2": float(integral), "reference_integrated_D2": exact_integral,
            "block_max_relative_L2_errors": block_errors,
            "checkpoints": "every proposed accepted interval, including t=0 power",
            "physical_power_quadrature": "trapezoidal", "final_state": q.tolist()}


def make_plan(reference, horizon, rho_relevant, input_fingerprints, policy=None):
    policy = policy or PlannerPolicy()
    required = ("operator", "equilibrium", "perturbation", "Krylov")
    if any(not input_fingerprints.get(k) for k in required):
        raise ValueError("All immutable input fingerprints required")
    exponential_error = reference.verify_modal_exponential(np.concatenate(([0.], np.geomspace(horizon*1e-8, horizon, 12))))
    if exponential_error > 1e-10:
        raise ValueError("Modal evaluation does not match direct reduced expm")
    base = candidate_blocks(reference, horizon, rho_relevant, policy)
    candidates = []
    for level in range(policy.max_candidate_refinements+1):
        blocks = refine_blocks(base, 2**level)
        if sum(b["steps"] for b in blocks) > policy.max_analysis_steps:
            break
        errors = verify_schedule(reference, blocks)
        passed = (errors["endpoint_relative_L2_error"] <= policy.endpoint_relative_tolerance and
                  errors["max_checkpoint_relative_L2_error"] <= policy.checkpoint_relative_tolerance and
                  errors["max_quadratic_energy_relative_error"] <= policy.energy_relative_tolerance and
                  errors["integrated_D2_relative_error"] <= policy.integrated_power_relative_tolerance)
        candidates.append({"uniform_refinement": level, "steps": sum(b["steps"] for b in blocks),
                           "passed": bool(passed), "predicted_errors": errors})
        if passed:
            for b, error in zip(blocks, errors["block_max_relative_L2_errors"]):
                b["predicted_max_error"] = error
                b["max_z_relevant"] = rho_relevant*b["dt"]
            dt0 = blocks[0]["dt"]
            z = rho_relevant*dt0
            result = {"algorithm_version": ALGORITHM_VERSION, "qualification_status": "linear_reduced_pass",
                "full_sparse_verification_required": True, "nonlinear_execution_authorized": False,
                "input_fingerprints": copy.deepcopy(input_fingerprints), "reduced_model_fingerprint": reference.fingerprint,
                "policy": asdict(policy), "policy_sha256": json_hash(asdict(policy)),
                "z_target": policy.z_target, "rho_relevant": rho_relevant, "tau_fast": 1/rho_relevant,
                "candidate_dt0": policy.z_target/rho_relevant, "dt0": dt0, "dt0_over_tau_fast": z,
                "max_BE_modal_amplification_error_first_block": float(abs(1/(1+z)-np.exp(-z))),
                "horizon": horizon, "blocks": blocks, "total_steps": sum(b["steps"] for b in blocks),
                "predicted_errors": errors, "candidate_history": candidates,
                "modal_vs_expm_relative_error": exponential_error}
            result["plan_sha256"] = json_hash(result)
            return result
    return {"qualification_status": "failed", "nonlinear_execution_authorized": False,
            "reason": "No candidate passed all predeclared linear gates within finite analysis guard",
            "candidate_history": candidates, "policy": asdict(policy), "input_fingerprints": input_fingerprints}


def validate_plan(plan, input_fingerprints):
    if plan.get("plan_sha256") != json_hash({k: v for k, v in plan.items() if k != "plan_sha256"}):
        raise ValueError("Plan integrity mismatch")
    if plan["input_fingerprints"] != input_fingerprints:
        raise ValueError("Plan operator/equilibrium/perturbation/Krylov fingerprint mismatch")
    if plan["policy_sha256"] != json_hash(plan["policy"]):
        raise ValueError("Plan policy fingerprint mismatch")


def verify_full_sparse_schedule(operator, reference, basis, blocks, monitor=None):
    """Check coupled sparse BE against the qualified lifted exponential.

    Computes physical D2 independently at every endpoint, including t=0.
    The analysis projection in operator.be_step removes roundoff mass only.
    No dense full operator/inverse and no energy-derived power integral.
    """
    basis = np.asarray(basis)
    if basis.shape != (operator.size, reference.size):
        raise ValueError("Incompatible Krylov lifting basis")
    q = basis @ reference.seed
    initial_norm = operator.norm(q)
    initial_energy = operator.energy(q)
    previous_power = operator.power(q)
    integral = 0.
    rows, block_rows = [], []
    try:
        for block in blocks:
            operator.clear_be_cache()
            local = []
            dt = block["dt"]
            for j in range(block["steps"]):
                t = block["t_start"]+(j+1)*dt
                q = operator.be_step(q, dt)
                exact = basis @ reference.state(t)
                power = operator.power(q)
                integral += dt*(previous_power+power)/2
                previous_power = power
                row = {"time": t, "dt": dt,
                    "relative_L2_error": operator.norm(q-exact)/max(operator.norm(exact), 1e-14*initial_norm),
                    "relative_energy_error": abs(operator.energy(q)-operator.energy(exact))/
                        max(abs(operator.energy(exact)), 1e-14*abs(initial_energy)),
                    "D2": power, "integrated_D2": integral,
                    "mass_domain": float(abs(operator.mass @ q)/operator.area)}
                rows.append(row); local.append(row)
            report = {"block": block["block"], "time": t, "dt": dt, "steps": block["steps"],
                "max_relative_L2_error": max(r["relative_L2_error"] for r in local),
                "max_relative_energy_error": max(r["relative_energy_error"] for r in local)}
            block_rows.append(report)
            if monitor:
                monitor(report)
    finally:
        operator.clear_be_cache()
    exact_integral = reference.integrated_power(blocks[-1]["t_end"])
    return {"endpoint_relative_L2_error": rows[-1]["relative_L2_error"],
        "max_checkpoint_relative_L2_error": max(r["relative_L2_error"] for r in rows),
        "max_quadratic_energy_relative_error": max(r["relative_energy_error"] for r in rows),
        "integrated_D2_relative_error": abs(integral/exact_integral-1),
        "linear_BE_integrated_D2": integral, "reference_integrated_D2": exact_integral,
        "max_mass_domain": max(r["mass_domain"] for r in rows),
        "blocks": block_rows, "history": rows}


class FrozenSchedule:
    """Production-compatible read-only schedule; no field-dependent API."""
    def __init__(self, plan, input_fingerprints, factor=1, analysis_only=False):
        validate_plan(plan, input_fingerprints)
        if not analysis_only and not plan.get("nonlinear_execution_authorized", False):
            raise ValueError("Unqualified plan cannot authorize nonlinear execution")
        self._blocks = tuple(tuple((k, v) for k, v in b.items()) for b in refine_blocks(plan["blocks"], factor))
        self.nsteps = sum(dict(b)["steps"] for b in self._blocks)
        self.plan_sha256 = plan["plan_sha256"]
        self.factor = factor

    @property
    def descriptor(self):
        return {"blocks": [dict(b) for b in self._blocks], "total_steps": self.nsteps,
                "plan_sha256": self.plan_sha256, "refinement_factor": self.factor,
                "switch_rule": "immutable fingerprinted all-BE blocks; no runtime adaptation"}

    def step(self, accepted_steps):
        if not 0 <= accepted_steps < self.nsteps:
            raise ValueError("No further planned interval")
        local = accepted_steps
        for data in self._blocks:
            block = dict(data)
            if local < block["steps"]:
                return {"dt": block["dt"], "time": block["t_start"]+(local+1)*block["dt"],
                        "phase": "spectral_plan_be", "scheme": "be"}
            local -= block["steps"]
        raise RuntimeError("Incomplete frozen plan")
