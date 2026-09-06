"""Exact cross-mesh quadrature: intersect both tensor-grid diagonals.

Every resulting integration triangle lies in ONE element of EACH mesh.
Degree-4 quadrature integrates (P2-P2)^2 and (DG1-DG1)^2 exactly.
No visualization grid, projection or smoothing is used for error norms.
"""

import numpy as np
from skfem.quadrature import get_quadrature_tri


def clip_halfplane(polygon, origin, normal):
    output = []
    for k, point in enumerate(polygon):
        previous = polygon[k-1]
        d, old = np.dot(point-origin, normal), np.dot(previous-origin, normal)
        inside, was_inside = d >= 0., old >= 0.
        if inside != was_inside:
            output.append(previous+old/(old-d)*(point-previous))
        if inside:
            output.append(point)
    return output


def diagonals(geometry):
    candidates = geometry.cell_triangles
    if np.any(candidates[:, 2] >= 0):
        raise ValueError("Overlay currently requires two triangles per tensor cell (even Nx)")
    result = []
    for first, second in candidates[:, :2]:
        shared = np.intersect1d(geometry.t[:, first], geometry.t[:, second])
        if len(shared) != 2:
            raise ValueError("Cell triangles must share one diagonal")
        a, b = geometry.p[:, shared].T
        edge = b-a
        result.append((a, np.array([-edge[1], edge[0]])))
    return result


def common_quadrature(a, b, order=4):
    xs = np.unique(np.r_[a.x, b.x, -.9, .9])
    zs = np.unique(np.r_[a.z, b.z, -1.])
    xs = xs[(xs >= a.x[0]) & (xs <= a.x[-1])]
    zs = zs[(zs >= a.z[0]) & (zs <= a.z[-1])]
    da, db = diagonals(a), diagonals(b)
    cx, cz = np.meshgrid((xs[:-1]+xs[1:])/2, (zs[:-1]+zs[1:])/2, indexing="ij")
    ca, cb = a.cell_indices(np.array([cx.ravel(), cz.ravel()])), b.cell_indices(np.array([cx.ravel(), cz.ravel()]))
    pieces = []
    n = 0
    for i in range(len(xs)-1):
        for j in range(len(zs)-1):
            polygons = [[np.array([xs[i], zs[j]]), np.array([xs[i+1], zs[j]]),
                         np.array([xs[i+1], zs[j+1]]), np.array([xs[i], zs[j+1]])]]
            for origin, normal in (da[ca[n]], db[cb[n]]):
                new = []
                for polygon in polygons:
                    for sign in (1., -1.):
                        cut = clip_halfplane(polygon, origin, sign*normal)
                        if len(cut) >= 3:
                            new.append(cut)
                polygons = new
            for polygon in polygons:
                for k in range(1, len(polygon)-1):
                    pieces.append([polygon[0], polygon[k], polygon[k+1]])
            n += 1
    vertices = np.asarray(pieces)
    determinants = abs(np.linalg.det((vertices[:, 1:]-vertices[:, :1]).transpose(0, 2, 1)))
    keep = determinants > 1e-25
    vertices, determinants = vertices[keep], determinants[keep]
    uv, weights = get_quadrature_tri(order)
    points = vertices[:, :1]+np.einsum("eij,jq->eqi", (vertices[:, 1:]-vertices[:, :1]).transpose(0, 2, 1), uv)
    weights = determinants[:, None]*weights
    return points.reshape(-1, 2).T, weights.ravel(), len(vertices)
