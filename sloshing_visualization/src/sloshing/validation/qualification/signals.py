"""Measured modal events and exact horizontal-line FEM integration."""

import numpy as np


def modal_events(times, signal, a=1., d=10., g=9.81):
    times, y = np.asarray(times), np.asarray(signal)
    crossings = []
    for i in range(len(times)-1):
        if y[i]*y[i+1] < 0:
            crossings.append(float(times[i]-y[i]*(times[i+1]-times[i])/(y[i+1]-y[i])))
    extrema = []
    for i in range(1, len(y)-1):
        if (y[i]-y[i-1])*(y[i+1]-y[i]) < 0:
            # Three-point local polynomial only to estimate an event time;
            # never used to change saved fields or particle trajectories.
            p = np.polyfit(times[i-1:i+2]-times[i], y[i-1:i+2], 2)
            offset = -p[1]/(2*p[0])
            t = times[i]+offset
            extrema.append({"time": float(t), "amplitude": float(np.polyval(p, offset)),
                            "kind": "maximum" if p[0] < 0 else "minimum"})
    opposite = next((r for r in extrema if r["amplitude"]*y[0] < 0), None)
    damping = []
    for kind in ("maximum", "minimum"):
        same = [e for e in extrema if e["kind"] == kind]
        for before, after in zip(same[:-1], same[1:]):
            ratio = abs(after["amplitude"]/before["amplitude"])
            damping.append({"from_time": before["time"], "to_time": after["time"], "amplitude_ratio": ratio,
                            "log_decrement": float(-np.log(ratio)),
                            "effective_decay_rate_per_s": float(-np.log(ratio)/(after["time"]-before["time"]))})
    k = np.pi/(2*a)
    return {"modal_basis": "sin(pi*x/(2*a)); diagnostic projection only",
            "inviscid_reference_period_s": float(2*np.pi/np.sqrt(g*k*np.tanh(k*d))),
            "zero_crossings_s": crossings, "first_zero_crossing_s": crossings[0] if crossings else None,
            "next_zero_crossing_s": crossings[1] if len(crossings) > 1 else None,
            "first_opposite_extremum": opposite, "extrema": extrema,
            "estimated_period_s": float(2*np.median(np.diff(crossings))) if len(crossings) > 1 else None,
            "half_cycle_period_estimates_s": list(2*np.diff(crossings)), "same_sign_damping": damping,
            "opposite_to_initial_amplitude_ratio": abs(opposite["amplitude"]/y[0]) if opposite else None}


def horizontal_quadrature(geometry, depth):
    # Split the horizontal line at every intersected FEM edge. P2 squared is
    # degree four on each remaining segment; three Gauss points are exact.
    edges = np.unique(np.sort(np.vstack((geometry.t[[0, 1]].T, geometry.t[[1, 2]].T,
                                        geometry.t[[0, 2]].T)), axis=1), axis=0)
    p0, p1 = geometry.p[:, edges[:, 0]].T, geometry.p[:, edges[:, 1]].T
    dz = p1[:, 1]-p0[:, 1]
    hit = (dz != 0) & (depth >= np.minimum(p0[:, 1], p1[:, 1])) & (depth <= np.maximum(p0[:, 1], p1[:, 1]))
    x = p0[hit, 0]+(depth-p0[hit, 1])/dz[hit]*(p1[hit, 0]-p0[hit, 0])
    x = np.unique(np.r_[geometry.x[0], x, geometry.x[-1]])
    nodes, weights = np.polynomial.legendre.leggauss(3)
    half = np.diff(x)/2
    queries = ((x[:-1]+x[1:])[:, None]/2+half[:, None]*nodes).ravel()
    return np.vstack((queries, np.full_like(queries, depth))), (half[:, None]*weights).ravel()
