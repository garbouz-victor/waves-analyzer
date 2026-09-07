"""Deterministic global-error greedy split/merge on a qualified reduced CH model.

No nonlinear state, runtime observations or soft acceptance penalties are inputs.
Local minimality is only with respect to the specified legal one-block moves.
"""
from dataclasses import dataclass
import copy
import time
import numpy as np
from .ch_timestep_planner import ReducedCHReference, verify_schedule, refine_blocks
from .ch_validation_policy import ValidationPolicy, LINEAR_LIMITS, linear_gates
from .linearized_ch import json_hash

VERSION = "global-hard-gate-blockwise-greedy-v1"


@dataclass(frozen=True)
class OptimizerGuard:
    cpu_s: float = 1800.
    evaluations: int = 5000
    rounds: int = 24
    max_splits: int = 3


def legal_blocks(blocks):
    previous_end = 0.
    for i, block in enumerate(blocks):
        if (block["steps"] <= 0 or block["steps"] != int(block["steps"]) or block["dt"] <= 0 or
                not np.isfinite([block["dt"], block["t_start"], block["t_end"]]).all()):
            return False
        scale = max(abs(block["t_end"]), 1e-300)
        if (abs(block["t_start"]-previous_end) > 1e-12*scale or
                abs(block["t_start"]+block["steps"]*block["dt"]-block["t_end"]) > 1e-12*scale):
            return False
        if i and block["dt"] > 2*blocks[i-1]["dt"]*(1+1e-12):
            return False
        previous_end = block["t_end"]
    return bool(blocks)


def original_from_uniform(plan):
    blocks = copy.deepcopy(plan["blocks"])
    if any(b["steps"] % 2 for b in blocks):
        raise ValueError("Uniform baseline is not an exact doubled original")
    for b in blocks:
        b["steps"] //= 2; b["dt"] *= 2
    if not legal_blocks(blocks):
        raise ValueError("Original blocks have invalid boundaries/growth")
    return blocks


def optimize(reference, original_blocks, input_fingerprints, *, reference_qualified,
             guard=OptimizerGuard(), candidate_callback=None):
    if not isinstance(reference, ReducedCHReference) or not reference_qualified:
        raise ValueError("Only a qualified reduced/reference linear model is accepted")
    if set(input_fingerprints) != {"operator", "equilibrium", "perturbation", "Krylov"}:
        raise ValueError("Require complete linear input fingerprints")
    if not legal_blocks(original_blocks):
        raise ValueError("Invalid original block schedule")
    started = time.process_time()
    history, cache = [], {}
    refinements = [0]*len(original_blocks)

    def blocks_at(levels):
        return [{**b, "dt": b["dt"]/2**level, "steps": b["steps"]*2**level}
                for b, level in zip(original_blocks, levels)]

    def evaluate(levels, action, block=None):
        if len(history) >= guard.evaluations or time.process_time()-started >= guard.cpu_s:
            raise TimeoutError("finite optimizer resource guard")
        key = tuple(levels)
        if key not in cache:
            cache[key] = verify_schedule(reference, blocks_at(levels))
        result = cache[key]
        errors = {k: float(result[k]) for k in LINEAR_LIMITS}
        passed = all(linear_gates(errors).values())
        if not np.isfinite(list(errors.values())).all():
            raise ValueError("Nonfinite reduced verification")
        row = {"candidate": len(history), "action": action, "block": block,
            "refinements": list(levels), "steps": sum(b["steps"] for b in blocks_at(levels)),
            "errors": errors, "passed": passed,
            "global_violation": max(0., max(errors[k]/v for k, v in LINEAR_LIMITS.items())-1.)}
        history.append(row)
        if candidate_callback:
            candidate_callback(copy.deepcopy(row))
        return row

    best_passing = None
    reason, local_minimal = "no passing candidate", False
    try:
        current = evaluate(refinements, "original")
        if current["passed"]:
            best_passing = current
        for _ in range(guard.rounds):
            if current["passed"]:
                break
            trials = []
            for i, level in enumerate(refinements):
                if level >= guard.max_splits:
                    continue
                trial_levels = refinements.copy(); trial_levels[i] += 1
                if not legal_blocks(blocks_at(trial_levels)):
                    continue
                trial = evaluate(trial_levels, "trial_split", i)
                if trial["passed"] and (best_passing is None or trial["steps"] < best_passing["steps"]):
                    best_passing = trial
                reduction = current["global_violation"]-trial["global_violation"]
                if reduction > 0:
                    trials.append((reduction/(trial["steps"]-current["steps"]), -i, trial))
            if not trials:
                reason = "no legal one-block split reduces global violation"
                break
            chosen = max(trials, key=lambda t: (t[0], t[1]))[2]
            refinements = chosen["refinements"].copy()
            current = evaluate(refinements, "accept_split", chosen["block"])
            if current["passed"]:
                best_passing = current
        if best_passing is not None:
            refinements = best_passing["refinements"].copy()
            while True:
                changed = False
                for i, level in enumerate(refinements):
                    if level == 0:
                        continue
                    trial_levels = refinements.copy(); trial_levels[i] -= 1
                    if not legal_blocks(blocks_at(trial_levels)):
                        continue
                    trial = evaluate(trial_levels, "trial_merge", i)
                    if trial["passed"]:
                        refinements = trial_levels
                        best_passing = evaluate(refinements, "accept_merge", i)
                        changed = True
                if not changed:
                    local_minimal = True
                    reason = "passing fixed point of legal one-block merges; not a global optimum"
                    break
    except TimeoutError as exc:
        reason = str(exc)+"; best-found passing schedule only, not proven optimal"
    policy = ValidationPolicy().record()
    record = {"optimizer_version": VERSION, "reference_fingerprint": reference.fingerprint,
        "input_fingerprints": dict(input_fingerprints), "policy": policy,
        "policy_sha256": json_hash(policy), "optimization_history": history,
        "qualification_status": "linear_reduced_pass" if best_passing else "failed",
        "reason": reason, "locally_minimal_under_legal_moves": local_minimal,
        "global_optimum_claimed": False,
        "guard": {"cpu_s": guard.cpu_s, "evaluations": guard.evaluations,
                  "rounds": guard.rounds, "max_splits": guard.max_splits}}
    if best_passing:
        blocks = blocks_at(best_passing["refinements"])
        record.update(blocks=blocks, total_steps=best_passing["steps"],
            predicted_errors=best_passing["errors"], dt0=blocks[0]["dt"],
            refinements=best_passing["refinements"], t_end=blocks[-1]["t_end"])
    record["plan_sha256"] = json_hash(record)
    return record
