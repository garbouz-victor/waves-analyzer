"""PW2 full-boundary engineering comparisons; never a graph eta(x) comparison."""
import json
import hashlib
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree
from skfem import Basis, ElementTriP2, ElementVector, MeshTri2

from .free_boundary_verify import (MODEL_ID, E0, SLOPE, _facet_quadrature,
                                   _hash_file, curve_samples, physical_contract)


def comparison_dependencies():
    directory = Path(__file__).parent
    return {name: _hash_file(directory/name) for name in
            ("free_boundary_compare.py", "free_boundary_verify.py")}


def _energy_observables(geometry, velocity, triangles, boundaries, intorder, grad_div_gamma=0.):
    mesh = MeshTri2(doflocs=geometry, t=triangles, _boundaries=boundaries)
    scalar = Basis(mesh, ElementTriP2(), intorder=intorder)
    vector = Basis(mesh, ElementVector(ElementTriP2()), intorder=intorder)
    field = vector.interpolate(velocity)
    kinetic = .5*float(np.sum(np.sum(np.asarray(field)**2, axis=0)*scalar.dx))
    potential = 9.81*(float(np.sum(np.asarray(scalar.global_coordinates())[1]*scalar.dx))+100.)
    strain_twice = field.grad+field.grad.swapaxes(0, 1)
    bulk = .005*float(np.sum(np.sum(strain_twice**2, axis=(0, 1))*scalar.dx))
    dilation = field.grad[0, 0]+field.grad[1, 1]
    divergence_squared = float(np.sum(dilation*dilation*scalar.dx))
    vertical = velocity[vector.split_indices()[1]]
    wall = 0.
    for q in _facet_quadrature(mesh, scalar, boundaries["left"], intorder):
        wall += .02*float(np.sum((vertical[q["dofs"]]@q["phi"])**2*q["measure"]))
    return [kinetic+potential, kinetic, potential, bulk, wall,
            grad_div_gamma*divergence_squared, float(np.sqrt(divergence_squared))]


def _load(path, heartbeat, until_time_s=None):
    with h5py.File(path, "r") as h:
        identity = json.loads(h.attrs["numerical_identity"])
        if identity.get("model_id") != MODEL_ID or identity.get("physical_contract") != physical_contract():
            raise ValueError("refinement_cannot_mix_linear_or_different_physics")
        topology = h["topology"]
        ids = topology["interface_nodes"][:]
        triangles = topology["triangles"][:]
        initial_geometry = topology["initial_geometry"][:]
        initial_mesh_digest = hashlib.sha256(initial_geometry.tobytes()+triangles.tobytes()).hexdigest()
        boundaries = {name: topology[name+"_facets"][:] for name in ("left", "right", "bottom", "surface")}
        time, surface, H, energies = [], [], [], []
        groups = [h["states"][key] for key in sorted(h["states"], key=int)]
        all_times = np.asarray([float(group.attrs["time_s"]) for group in groups])
        if (len(all_times) < 2 or not np.isfinite(all_times).all() or all_times[0] != 0.
                or np.min(np.diff(all_times)) <= 0.):
            raise ValueError("comparison_requires_ordered_accepted_native_history")
        full_native_horizon = float(all_times[-1]); full_native_states = len(groups)
        if until_time_s is not None:
            matches = np.flatnonzero(abs(all_times-until_time_s) <= 1e-12)
            if len(matches) != 1 or matches[0] == 0:
                raise ValueError("explicit_prefix_end_must_match_one_accepted_state_in_each_native")
            groups = groups[:int(matches[0])+1]
        for index, group in enumerate(groups):
            if not bool(group.attrs.get("accepted", False)):
                raise ValueError("comparison_prefix_contains_nonaccepted_native_state")
            geometry = group["geometry"][:]
            time.append(float(group.attrs["time_s"]))
            surface.append(geometry[:, ids])
            H.append(json.loads(group.attrs["coating"])["H"])
            energies.append(_energy_observables(geometry, group["velocity"][:], triangles, boundaries,
                                               identity["controls"]["intorder"],
                                               identity["controls"].get("grad_div_gamma_m2_s", 0.)))
            if index % 20 == 0:
                heartbeat()
        markers = json.loads(groups[-1].attrs["coating"])["markers"]
        for marker in markers:
            marker["source_time_s"] = time[int(marker["state_id"].rsplit(":", 1)[1])]
        if len(time) < 2 or time[0] != 0. or np.min(np.diff(time)) <= 0.:
            raise ValueError("comparison_requires_ordered_accepted_native_history")
        if np.max(np.diff(time)) > identity["controls"]["dt"]+1e-12:
            raise ValueError("native_step_exceeds_declared_dt")
    return {"identity": identity, "times": np.asarray(time), "surface": np.asarray(surface),
            "H": np.asarray(H), "energy": np.asarray(energies), "markers": markers,
            "initial_mesh_digest": initial_mesh_digest, "actual_triangles": triangles.shape[1],
            "full_native_horizon_s": full_native_horizon, "full_native_states": full_native_states}


def _at(data, time_s):
    # A comparison endpoint may differ from its committed native clock by
    # <=1e-12 s. Use that actual state, never extrapolate an invented endpoint.
    if abs(time_s-data["times"][-1]) <= 1e-12:
        return [data[name][-1] for name in ("surface", "H", "energy")]
    right = min(np.searchsorted(data["times"], time_s, side="right"), len(data["times"])-1)
    left = max(right-1, 0)
    fraction = 0. if left == right else (time_s-data["times"][left])/(data["times"][right]-data["times"][left])
    return [(1-fraction)*data[name][left]+fraction*data[name][right] for name in ("surface", "H", "energy")]


def _dense(line, maximum_spacing):
    pieces = []; actual_spacing = 0.
    for start, end in zip(line[:-1], line[1:]):
        length = float(np.linalg.norm(end-start))
        count = max(1, int(np.ceil(length/maximum_spacing)))
        pieces.append(start[None, :]+(np.arange(count)/count)[:, None]*(end-start)[None, :])
        actual_spacing = max(actual_spacing, length/count)
    return np.vstack(pieces+[line[-1:]]), actual_spacing/2.


def _directed(points, target):
    """Original sampled point-to-segment distance with conservative pruning.

    A nearest target vertex gives an upper bound U on the segment distance.
    If a segment with midpoint m and half-length r can improve that bound,
    |q-m| <= U+r.  The exact KD radius U+max(r), padded outwards for roundoff,
    therefore retains every possible nearest segment.  No nearest-midpoint
    approximation is used.  The actual distance arithmetic remains unchanged.
    """
    start = target[:-1]
    edge = target[1:]-start
    denominator = np.maximum(np.sum(edge*edge, axis=1), 1e-300)
    def brute(chunk):
        point = chunk[:, None, :]
        fraction = np.clip(np.sum((point-start)*edge, axis=2)/denominator, 0., 1.)
        difference = point-start-fraction[:, :, None]*edge
        return float(np.sqrt(np.max(np.min(np.sum(difference**2, axis=2), axis=1))))
    def brute_all():
        maximum = 0.
        for first in range(0, len(points), 128):
            maximum = max(maximum, brute(points[first:first+128]))
        return maximum
    if len(points) == 0:
        return 0.
    if (not len(edge) or not np.isfinite(points).all() or not np.isfinite(target).all()
            or not np.isfinite(edge).all() or not np.isfinite(denominator).all()):
        return brute_all()
    # Half sums avoid an unnecessary overflowing endpoint addition.  Very
    # small coordinates may receive generous padding; that only adds work.
    midpoint = .5*start+.5*target[1:]
    half_length = .5*np.sqrt(np.sum(edge*edge, axis=1))
    maximum_radius = float(np.max(half_length))
    coordinate_scale = max(1., float(np.max(abs(target))), float(np.max(abs(points))))
    try:
        vertices = cKDTree(target)
        segments = cKDTree(midpoint)
        maximum = 0.
        for first in range(0, len(points), 128):
            chunk = points[first:first+128]
            upper = vertices.query(chunk, k=1, eps=0., p=2)[0]
            radius = upper+maximum_radius
            padding = 64*np.finfo(float).eps*(coordinate_scale+radius)
            radius = np.nextafter(radius+padding, np.inf)
            if not np.isfinite(radius).all():
                maximum = max(maximum, brute(chunk))
                continue
            candidates = segments.query_ball_point(chunk, radius, p=2, eps=0., return_sorted=True)
            counts = np.asarray([len(indices) for indices in candidates], dtype=int)
            if np.any(counts == 0):
                maximum = max(maximum, brute(chunk))
                continue
            selected = np.concatenate(candidates).astype(np.intp, copy=False)
            queries = np.repeat(chunk, counts, axis=0)
            fraction = np.clip(np.sum((queries-start[selected])*edge[selected], axis=1)
                               /denominator[selected], 0., 1.)
            difference = queries-start[selected]-fraction[:, None]*edge[selected]
            squared = np.sum(difference*difference, axis=1)
            if not np.isfinite(squared).all():
                maximum = max(maximum, brute(chunk))
                continue
            offsets = np.r_[0, np.cumsum(counts[:-1])]
            nearest = np.minimum.reduceat(squared, offsets)
            maximum = max(maximum, float(np.sqrt(np.max(nearest))))
        return maximum
    except (ValueError, OverflowError, FloatingPointError):
        # No approximate result is returned if tree arithmetic is unsuitable.
        return brute_all()


def curve_distance(first, second, *, tessellation_m=1e-6, maximum_spacing_m=.001):
    """Two-sided Hausdorff bracket with explicit spatial sampling error.

    Distance to any closed set is 1-Lipschitz. Adding half the source sample
    spacing and both P2 chord bounds controls maxima between source samples.
    """
    one, error1 = curve_samples(first, np.arange(first.shape[1]), tessellation_m)
    two, error2 = curve_samples(second, np.arange(second.shape[1]), tessellation_m)
    points1, gap1 = _dense(one, maximum_spacing_m)
    points2, gap2 = _dense(two, maximum_spacing_m)
    forward, backward = _directed(points1, two), _directed(points2, one)
    sampled = max(forward, backward)
    return {"sampled_m": sampled, "lower_bound_m": max(0., sampled-error1-error2),
            "upper_bound_m": max(forward+gap1, backward+gap2)+error1+error2,
            "summed_tessellation_bound_m": error1+error2,
            "sampling_half_gap_bound_m": max(gap1, gap2)}


def _marker_differences(first, second):
    """Compare identified accepted peaks without creating absent records.

    IDs encode chronological resolved records, not a fitted correspondence
    between interior crests.  A changed record catalog remains explicit.
    """
    def records(markers):
        previous = None
        result = {}
        for marker in markers:
            if marker["side"] != "L":
                continue
            if previous is not None:
                result[marker["marker_id"]] = {
                    **marker, "record_excess_m": marker["height_m"]-previous}
            previous = marker["height_m"]
        return result
    left, right = records(first), records(second)
    matched = []
    for name in left:
        if name not in right:
            continue
        one, two = left[name], right[name]
        difference = abs(one["height_m"]-two["height_m"])
        matched.append({"marker_id": name,
                        "height_difference_m": difference,
                        "source_time_difference_s": abs(one["source_time_s"]-two["source_time_s"]),
                        "confirmation_time_difference_s": abs(one["created_at_s"]-two["created_at_s"]),
                        "record_excess_m": [one["record_excess_m"], two["record_excess_m"]],
                        "record_excess_not_larger_than_height_difference":
                            min(one["record_excess_m"], two["record_excess_m"]) <= difference})
    return {"matched_chronological_records": matched,
            "only_in_first": [name for name in left if name not in right],
            "only_in_second": [name for name in right if name not in left],
            "detector_floor_m": 1e-6,
            "qualification": "The detector floor is not an uncertainty bound. Differences are measured between levels; a small record excess or changed catalog is explicitly not robust record identity evidence."}


def compare_runs(first, second, *, kind=None, require_target=True, sample_dt=None,
                 maximum_spacing_m=.001, heartbeat=lambda: None, until_time_s=None):
    """Compare every union accepted time by default; optional explicit cadence.

    The parameterized native P2 curves are interpolated in TIME, not sorted in
    x. Remeshing in this method preserves all surface nodes. Bounds here are
    spatial geometric bounds at comparison times, not rigorous temporal error.
    An explicit until_time_s selects an accepted common prefix read-only. Its
    endpoint must exist in BOTH natives (within clock roundoff); no interpolated
    terminal state or future marker is substituted. Partial prefixes cannot
    qualify a release. Defaults remain full-native, 0-5 s comparisons.
    """
    first, second = Path(first), Path(second)
    if until_time_s is not None:
        try:
            until_time_s = float(until_time_s)
        except (TypeError, ValueError) as error:
            raise ValueError("explicit_prefix_end_must_be_positive_finite_seconds") from error
        if not np.isfinite(until_time_s) or until_time_s <= 0.:
            raise ValueError("explicit_prefix_end_must_be_positive_finite_seconds")
        if require_target and abs(until_time_s-5.) > 1e-12:
            raise ValueError("explicit_partial_prefix_cannot_qualify_0_to_5s")
        if sample_dt is not None:
            raise ValueError("explicit_prefix_requires_union_of_all_accepted_times_not_sample_dt")
    one, two = _load(first, heartbeat, until_time_s), _load(second, heartbeat, until_time_s)
    if until_time_s is None and abs(one["times"][-1]-two["times"][-1]) > 1e-12:
        raise ValueError("refinement_horizons_differ")
    if one["identity"]["source_hashes"] != two["identity"]["source_hashes"]:
        raise ValueError("different_numerical_sources_are_not_a_controlled_dt_h_comparison")
    c1, c2 = one["identity"]["controls"], two["identity"]["controls"]
    formulation = ("hydrostatic_split", "skew_divergence", "rezone_interval_s", "rezone_strategy", "intorder",
                   "grad_div_gamma_m2_s")
    if any(c1.get(key, False if key.endswith("split") or key.endswith("divergence") else 0.)
           != c2.get(key, False if key.endswith("split") or key.endswith("divergence") else 0.) for key in formulation):
        raise ValueError("different_numerical_formulations_are_not_dt_h_refinement")
    mesh_keys = ("nx", "nz", "x_grading", "z_grading", "local_right_levels")
    if kind == "temporal":
        if any(c1[key] != c2[key] for key in mesh_keys) or not c2["dt"] < c1["dt"] or one["initial_mesh_digest"] != two["initial_mesh_digest"]:
            raise ValueError("temporal_refinement_requires_identical_mesh_and_smaller_dt")
    elif kind == "spatial":
        if c1["dt"] != c2["dt"] or two["actual_triangles"] <= one["actual_triangles"] or not (c2["nx"] > c1["nx"] or c2["nz"] > c1["nz"]
                                       or c2["local_right_levels"] > c1["local_right_levels"]):
            raise ValueError("spatial_refinement_requires_same_dt_and_finer_mesh")
    elif kind is not None:
        raise ValueError("comparison_kind_must_be_temporal_or_spatial")
    end = float(one["times"][-1]) if until_time_s is None else until_time_s
    excludes_native_future = any(len(data["times"]) < data["full_native_states"] for data in (one, two))
    if require_target and excludes_native_future:
        raise ValueError("explicit_partial_prefix_cannot_qualify_0_to_5s")
    if require_target and abs(end-5.) > 1e-12:
        raise ValueError("refinement_is_not_0_to_5s")
    if sample_dt is None:
        samples = np.unique(np.round(np.r_[one["times"], two["times"]], 12))
        if until_time_s is not None:
            # Only align the two already-matched native endpoint clocks.
            samples[abs(samples-end) <= 1e-12] = end
            samples = np.unique(samples)
    else:
        if not sample_dt > 0.:
            raise ValueError("sample_dt_must_be_positive")
        samples = np.unique(np.r_[np.arange(0., end, sample_dt), end])
    maxima = {"surface_sampled_m": 0., "surface_upper_bound_m": 0.,
              "left_contact_m": 0., "left_H_m": 0., "P2_drift_m": 0.}
    witnesses = {}; energy_difference = np.zeros(7)
    maximum_tessellation = 0.
    for index, time_s in enumerate(samples):
        X1, H1, E1 = _at(one, time_s); X2, H2, E2 = _at(two, time_s)
        distance = curve_distance(X1, X2, maximum_spacing_m=maximum_spacing_m)
        values = {"surface_sampled_m": distance["sampled_m"], "surface_upper_bound_m": distance["upper_bound_m"],
                  "left_contact_m": abs(X1[1, 0]-X2[1, 0]), "left_H_m": abs(H1[0]-H2[0]),
                  "P2_drift_m": max(abs(X1[1, -1]-SLOPE), abs(X2[1, -1]-SLOPE))}
        for name, value in values.items():
            if value >= maxima[name]:
                maxima[name] = float(value); witnesses[name] = float(time_s)
        maximum_tessellation = max(maximum_tessellation, distance["summed_tessellation_bound_m"])
        energy_difference = np.maximum(energy_difference, abs(E1-E2))
        if index % 10 == 0:
            heartbeat()
    span = 2*SLOPE
    energy_names = ("total_energy", "kinetic_energy", "potential_energy", "bulk_dissipation_rate", "left_wall_dissipation_rate",
                    "numerical_grad_div_dissipation_rate")
    normalized = {name: value/span for name, value in maxima.items() if name != "P2_drift_m"}
    return {"kind": kind, "passed": max(maxima[name] for name in ("surface_upper_bound_m", "left_contact_m", "left_H_m"))/span <= .05
                         and maxima["P2_drift_m"] <= 1e-12,
        "first_native": str(first), "second_native": str(second),
        "first_sha256": _hash_file(first), "second_sha256": _hash_file(second),
        "comparison_code_sha256": _hash_file(__file__), "time_end_s": end,
        "require_target": bool(require_target),
        "explicit_common_interval": until_time_s is not None,
        "requested_until_time_s": until_time_s,
        "prefix_excludes_native_future_states": excludes_native_future,
        "comparison_scope": "EXPLICIT_ACCEPTED_COMMON_PREFIX_DIAGNOSTIC" if until_time_s is not None else "FULL_NATIVE_TRAJECTORIES",
        "full_native_horizons_s": [one["full_native_horizon_s"], two["full_native_horizon_s"]],
        "full_native_state_counts": [one["full_native_states"], two["full_native_states"]],
        "compared_native_end_times_s": [float(one["times"][-1]), float(two["times"][-1])],
        "compared_native_state_counts": [len(one["times"]), len(two["times"])],
        "full_5s_interval_compared": bool(not excludes_native_future and abs(end-5.) <= 1e-12),
        "prefix_endpoint_policy": "actual accepted state in each native within 1e-12 s of requested clock; exact native endpoint fields, no synthetic terminal state" if until_time_s is not None else None,
        "comparison_dependencies": comparison_dependencies(),
        "actual_triangles": [one["actual_triangles"], two["actual_triangles"]],
        "comparison_times": len(samples), "common_time_policy": (
            "union_of_all_accepted_times_in_explicit_prefix" if until_time_s is not None else
            "union_of_all_accepted_times" if sample_dt is None else f"explicit_{sample_dt:g}s_cadence"),
        "maximum_native_step_s": [float(np.max(np.diff(one["times"]))), float(np.max(np.diff(two["times"])))],
        "maxima": maxima, "normalized_by_initial_height_span": normalized, "maximum_at_time_s": witnesses,
        "summed_P2_tessellation_bound_m": maximum_tessellation,
        "distance_spatial_bound": "two directed point-to-segment distances plus half source sample gaps and both exact P2 chord-error bounds",
        "temporal_interpolation": "linear interpolation of native material surface coefficients; no rigorous temporal interpolation error bound claimed",
        "energy_observable_differences_per_density": dict(zip(energy_names, energy_difference.tolist())),
        "accepted_strong_divergence_L2_difference": float(energy_difference[6]),
        "energy_differences_over_E0": dict(zip(energy_names[:3], (energy_difference[:3]/E0).tolist())),
        "first_markers": one["markers"], "second_markers": two["markers"],
        "peak_differences": _marker_differences(one["markers"], two["markers"]),
        "record_identity_agreement": [m["marker_id"] for m in one["markers"]] == [m["marker_id"] for m in two["markers"]],
        "formal_convergence_order_claimed": False, "physical_error_bound_claimed": False}


comparison_metrics = compare_runs
