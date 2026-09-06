"""A new interpretation layer; historical STEP 1.6 decisions are immutable."""

import json
from pathlib import Path
import numpy as np

CONTACT_WARNING = "Pinned-contact region: not physically resolved by this model"
CORNER_WARNING = "Unresolved contact corner; solid walls below it obey no-slip"
FIXED_DOMAIN = "Velocity: fixed reference domain z <= 0. Surface: first-order interface displacement."
LIMITATIONS = [CONTACT_WARNING, FIXED_DOMAIN,
    "No moving contact line, wetting/dewetting, film deposition or wall drainage",
    "Pointwise corner omega and slope are not converged; fine is not an exact continuum solution",
    "No validated nonlinear particle drift or mixing prediction",
    "Wall-vorticity field difference is about 5.7%; bulk velocity difference about 0.0793%",
    "Taylor-Hood is weakly, not pointwise, divergence-free",
    "Display grids and color statistics are not scientific error norms"]


def load_scope(path=None):
    path = Path(path) if path else Path(__file__).resolve().parents[3]/"configs/animation_scope.json"
    scope = json.loads(path.read_text())
    if not 0 < scope["bulk_abs_x_max_m"] < scope["transition_abs_x_max_m"] < scope["domain_abs_x_max_m"]:
        raise ValueError("Invalid interpretation zones")
    return scope


def region_codes(x, scope):
    """0 bulk, 1 transition, 2 contact; boundary belongs to inner region."""
    x = np.abs(np.asarray(x))
    return np.where(x > scope["transition_abs_x_max_m"], 2,
                    np.where(x > scope["bulk_abs_x_max_m"], 1, 0))


def magnify_positions(initial, actual, factor):
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("Magnification must be a fixed positive scalar")
    return np.asarray(initial) + factor*(np.asarray(actual)-np.asarray(initial))


def magnify_surface(eta, factor=1.):
    return np.asarray(eta)*factor


def compatible_cache(stored, expected):
    if stored != expected:
        raise RuntimeError("Display cache/source fingerprint or scope changed; use a new --output directory")


def event_frames(times, modal):
    events = [("t0", "Release", 0.),
              ("first_zero", "First modal crossing", modal["first_zero_crossing_s"]),
              ("opposite_extremum", "Opposite modal extremum", modal["first_opposite_extremum"]["time"]),
              ("next_zero", "Next modal crossing", modal["next_zero_crossing_s"])]
    same = next(e for e in modal["extrema"] if e["kind"] == "maximum")
    events.append(("next_positive", "Next positive modal extremum", same["time"]))
    late = next(e for e in modal["extrema"] if e["time"] > 3.8)
    events.append(("late_cycle", "Damped late extremum", late["time"]))
    return [{"name": n, "phase": label, "event_time_s": t,
             "index": int(np.argmin(abs(np.asarray(times)-t))),
             "selected_time_s": float(times[np.argmin(abs(np.asarray(times)-t))])}
            for n, label, t in events]


def phase_label(t, modal, interval=.035):
    for z in modal["zero_crossings_s"]:
        if abs(t-z) <= interval:
            return "Modal crossing: compare kinetic energy"
    for e in modal["extrema"]:
        if abs(t-e["time"]) <= interval:
            return "Modal extremum / overshoot"
    if t < modal["first_zero_crossing_s"]:
        return "Release / acceleration" if t < .35*modal["first_zero_crossing_s"] else "Approaching first crossing"
    if t < modal["first_opposite_extremum"]["time"]:
        return "Inertia carries the surface beyond equilibrium"
    if t < modal["next_zero_crossing_s"]:
        return "Return toward equilibrium: watch flow reversal"
    return "Damped cycle: fixed scales preserve amplitude loss"
