"""Extrema of the actual P2 polynomial and certified triangle normal widths.

Reference nodes: (0,0), (1,0), (0,1), (.5,0), (.5,.5), (0,.5).
No raster or angular sampling participates in either certificate. Floating
point bounds are rounded outward; ill-conditioned cases explicitly fall back
to Bernstein (range) or diameter (width), never to a sampled hard gate.
"""
import numpy as np

ALGORITHM_VERSION = "triangle-p2-extrema-gradient-cone-v2"
_EPS = np.finfo(float).eps
_VERTICES = np.array([[0., 0.], [1., 0.], [0., 1.]])


def coefficients(values):
    v = np.asarray(values, dtype=np.longdouble)
    if v.shape[-1] != 6 or not np.isfinite(v).all():
        raise ValueError("Six finite reference P2 nodal values required")
    f = v[..., 0]
    a = 2*(v[..., 1]+f-2*v[..., 3])
    c = 2*(v[..., 2]+f-2*v[..., 5])
    d, e = v[..., 1]-f-a, v[..., 2]-f-c
    b = 4*(v[..., 4]-v[..., 3]-v[..., 5]+f)
    return np.stack((a, b, c, d, e, f), axis=-1)


def evaluate(coefs, points):
    a, b, c, d, e, f = np.moveaxis(np.asarray(coefs), -1, 0)
    x, y = np.moveaxis(np.asarray(points), -1, 0)
    return ((a*x+b*y+d)*x+(c*y+e)*y+f)


def physical_gradients(triangles, values, points, degree=2):
    """Affine gradients from actual nodal coefficients in extended precision.

    This avoids evaluating a sum of large basis gradients for a nearly constant
    field. In particular an exactly constant coefficient field has exactly zero
    gradient before the width fallback decision.
    """
    tri=np.asarray(triangles,dtype=np.longdouble)
    points=np.asarray(points,dtype=np.longdouble)
    v=np.asarray(values,dtype=np.longdouble)
    if degree==2:
        a,b,c,d,e,_=coefficients(v).T
        gx=2*a[:,None]*points[None,:,0]+b[:,None]*points[None,:,1]+d[:,None]
        gy=b[:,None]*points[None,:,0]+2*c[:,None]*points[None,:,1]+e[:,None]
    elif degree==1:
        gx=np.broadcast_to((v[:,1]-v[:,0])[:,None],(len(v),len(points)))
        gy=np.broadcast_to((v[:,2]-v[:,0])[:,None],(len(v),len(points)))
    else:
        raise ValueError("Only CG1/CG2 triangle gradients supported")
    first,second=tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]
    det=first[:,0]*second[:,1]-first[:,1]*second[:,0]
    if (det==0).any():
        raise ValueError("Degenerate triangle mapping")
    return np.asarray(np.stack(((gx*second[:,1,None]-gy*first[:,1,None])/det[:,None],
                               (-gx*second[:,0,None]+gy*first[:,0,None])/det[:,None]),axis=-1),float)


def polynomial_ranges(values, degree=2):
    """Batch actual extrema with witnesses; fallback bounds have no witness.

    Returned min/max include a 128-eps coefficient-scale outward roundoff
    allowance. `min_value`/`max_value` are unpadded candidate values. A nearly
    singular but nonzero Hessian determinant is not treated as exactly zero.
    A truly singular Hessian needs only boundary extrema (stationary lines
    reach the boundary). All indices refer to reference coordinates.
    """
    v = np.asarray(values, dtype=float)
    if v.ndim != 2 or degree not in (1, 2) or v.shape[1] != (3 if degree == 1 else 6):
        raise ValueError("Expected (cells, 3) CG1 or (cells, 6) CG2 values")
    if not np.isfinite(v).all():
        raise ValueError("Nonfinite polynomial")
    n = len(v)
    points = np.broadcast_to(_VERTICES, (n, 3, 2)).copy()
    candidates = v[:, :3].astype(np.longdouble)
    types = ["vertex"]*3
    fallback = np.zeros(n, dtype=bool)
    if degree == 2:
        co = coefficients(v)
        a, b, c, d, e, f = co.T
        scale = np.maximum(1., np.abs(co).sum(axis=1))
        edge_A = (a, a-b+c, c)
        edge_B = (d, -2*a+b-d+e, e)
        starts = (_VERTICES[0], _VERTICES[1], _VERTICES[0])
        directions = (_VERTICES[1], np.array([-1., 1.]), _VERTICES[2])
        for A, B, start, direction in zip(edge_A, edge_B, starts, directions):
            t = np.zeros(n, dtype=np.longdouble)
            np.divide(-B, 2*A, out=t, where=np.abs(A) > 64*_EPS*scale)
            valid = (np.abs(A) > 64*_EPS*scale) & (t > 0) & (t < 1)
            # Invalid candidates evaluate at a vertex; no out-of-domain witness.
            p = start+np.where(valid, t, 0.)[:, None]*direction
            points = np.concatenate((points, p[:, None, :]), axis=1)
            candidates = np.column_stack((candidates, evaluate(co, p)))
            types.append("edge")
        det = 4*a*c-b*b
        hscale = np.maximum(np.maximum(abs(2*a), abs(2*c)), abs(b))
        singular = det == 0
        uncertain = (~singular) & (abs(det) <= 1e-12*hscale*hscale)
        fallback |= uncertain
        good = ~(singular | uncertain)
        x, y = np.zeros(n, dtype=np.longdouble), np.zeros(n, dtype=np.longdouble)
        np.divide(b*e-2*c*d, det, out=x, where=good)
        np.divide(b*d-2*a*e, det, out=y, where=good)
        valid = good & (x >= 0) & (y >= 0) & (x+y <= 1)
        p = np.column_stack((np.where(valid, x, 0.), np.where(valid, y, 0.)))
        points = np.concatenate((points, p[:, None, :]), axis=1)
        candidates = np.column_stack((candidates, evaluate(co, p)))
        types.append("interior")
        edge = 2*v[:, 3:6]-.5*v[:, [0, 1, 2]]-.5*v[:, [1, 2, 0]]
        bernstein = np.column_stack((v[:, :3], edge))
    else:
        scale = np.maximum(1., abs(v).sum(axis=1))
    imin, imax = candidates.argmin(axis=1), candidates.argmax(axis=1)
    rows = np.arange(n)
    lower, upper = candidates[rows, imin].astype(float), candidates[rows, imax].astype(float)
    min_point, max_point = points[rows, imin].astype(float), points[rows, imax].astype(float)
    min_type, max_type = np.array(types)[imin], np.array(types)[imax]
    min_value, max_value = lower.copy(), upper.copy()
    if degree == 2:
        lower = np.where(fallback, bernstein.min(axis=1), lower)
        upper = np.where(fallback, bernstein.max(axis=1), upper)
    allowance = np.asarray(128*_EPS*scale, dtype=float) if degree == 2 else np.zeros(n)
    return {"min": np.nextafter(lower-allowance, -np.inf) if degree == 2 else lower,
            "max": np.nextafter(upper+allowance, np.inf) if degree == 2 else upper,
            "min_value": min_value, "max_value": max_value,
            "min_point": min_point, "max_point": max_point,
            "min_type": np.where(fallback, "fallback", min_type),
            "max_type": np.where(fallback, "fallback", max_type),
            "exact": ~fallback,
            "status": np.where(fallback, "certified_conservative_fallback", "exact"),
            "fallback_reason": np.where(fallback, "nearly_singular_hessian", ""),
            "roundoff_allowance": allowance}


def certified_normal_widths(triangles, vertex_gradients):
    """Maximum projected width over the complete affine-gradient angular cone.

    Nonzero affine gradients lie in the convex hull of the three vertex
    gradients. Outside the origin this hull subtends an arc shorter than pi.
    For every oriented vertex pair, maximize d.n over that arc: endpoints or
    n parallel to d. A near-origin hull uses the triangle diameter instead.
    Widths are enlarged by a scale-relative floating point safety allowance.
    """
    tri, g = np.asarray(triangles, float), np.asarray(vertex_gradients, float)
    if tri.shape != g.shape or tri.ndim != 3 or tri.shape[1:] != (3, 2):
        raise ValueError("Triangles and vertex gradients must have shape (cells,3,2)")
    if not np.isfinite(tri).all() or not np.isfinite(g).all():
        raise ValueError("Nonfinite geometry or gradient")
    edges = tri-np.roll(tri, 1, axis=1)
    diameter = np.linalg.norm(edges, axis=-1).max(axis=1)
    if (diameter <= 0).any():
        raise ValueError("Degenerate triangle")
    norms = np.linalg.norm(g, axis=-1)
    scale = norms.max(axis=1)
    near = norms.min(axis=1) <= 1e-10*np.maximum(scale, 1e-300)
    angles = np.sort(np.arctan2(g[..., 1], g[..., 0]), axis=1)
    gaps = np.diff(np.column_stack((angles, angles[:, 0]+2*np.pi)), axis=1)
    gap_index = gaps.argmax(axis=1)
    rows = np.arange(len(tri))
    span = 2*np.pi-gaps[rows, gap_index]
    start = angles[rows, (gap_index+1) % 3]
    # Distance to hull edges also detects a narrow, ill-conditioned near miss.
    for i in range(3):
        u, v = g[:, i], g[:, (i+1) % 3]
        d = v-u
        denom = (d*d).sum(axis=1)
        t = np.zeros(len(tri))
        np.divide(-(u*d).sum(axis=1), denom, out=t, where=denom > 0)
        t = np.minimum(1., np.maximum(0., t))
        distance = np.linalg.norm(u+t[:, None]*d, axis=1)
        near |= distance <= 1e-10*np.maximum(scale, 1e-300)
    inside = span >= np.pi-1e-10
    fallback = near | inside
    # Widen cone endpoints to account for atan2/arc roundoff.
    start -= 1e-12
    span += 2e-12
    n0 = np.column_stack((np.cos(start), np.sin(start)))
    n1 = np.column_stack((np.cos(start+span), np.sin(start+span)))
    width = np.zeros(len(tri))
    for d in (edges[:, 0], edges[:, 1], edges[:, 2],
              -edges[:, 0], -edges[:, 1], -edges[:, 2]):
        relative_angle = (np.arctan2(d[:, 1], d[:, 0])-start) % (2*np.pi)
        candidate = np.maximum((d*n0).sum(axis=1), (d*n1).sum(axis=1))
        candidate = np.where(relative_angle <= span, np.linalg.norm(d, axis=1), candidate)
        width = np.maximum(width, candidate)
    width = np.minimum(width, diameter)
    width = np.where(fallback, diameter, width)
    width = np.nextafter(width+64*_EPS*diameter, np.inf)
    return {"width": width, "diameter": diameter, "diameter_fallback": fallback,
            "fallback_reason": np.where(inside, "gradient_hull_contains_origin",
                                   np.where(near, "gradient_hull_near_origin", "")),
            "cone_start_rad": start, "cone_span_rad": span}
