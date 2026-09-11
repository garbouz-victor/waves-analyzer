"""Actual right-wall-facing branch gap versus adjacent element normal width.

Diagnostic only: the gap vanishes at the material terminating contact P2. A
minimum gap/cell-width ratio is not by itself a continuum/model or release FAIL.
No molecular thickness, prescribed film shape, or cropped branch is introduced.
"""
import numpy as np
from skfem import Basis, ElementTriP2, MeshTri2

from .free_boundary_verify import SLOPE, curve_samples


def _quadratic_range(nodal):
    points = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
    polynomial = np.array([[1., x, y, x*x, x*y, y*y] for x, y in points])
    constant, x, y, xx, xy, yy = np.linalg.solve(polynomial, nodal)
    values = [float(value) for value in nodal[:3]]
    hessian = np.array([[2*xx, xy], [xy, 2*yy]])
    linear = np.array([x, y])
    candidates = []
    for start, direction in ((np.array([0., 0.]), np.array([1., 0.])),
                             (np.array([0., 0.]), np.array([0., 1.])),
                             (np.array([1., 0.]), np.array([-1., 1.]))):
        a = .5*direction@hessian@direction
        b = direction@(hessian@start+linear)
        if a != 0.:
            q = -b/(2*a)
            if 0. < q < 1.:
                candidates.append(start+q*direction)
    if abs(np.linalg.det(hessian)) > 1e-28*max(np.linalg.norm(hessian)**2, 1e-30):
        stationary = np.linalg.solve(hessian, -linear)
        if np.min(stationary) > 0. and np.sum(stationary) < 1.:
            candidates.append(stationary)
    for point in candidates:
        values.append(float(constant+linear@point+.5*point@hessian@point))
    return min(values), max(values)


def right_branch_gap(geometry, triangles, boundaries, interface_ids, *, tessellation_m=1e-6):
    mesh = MeshTri2(doflocs=geometry, t=triangles, _boundaries=boundaries)
    scalar = Basis(mesh, ElementTriP2())
    facets = {tuple(sorted(mesh.facets[:, facet])): int(facet) for facet in boundaries["surface"]}
    pieces, owner = [], []
    for edge in range((len(interface_ids)-1)//2):
        nodes = geometry[:, interface_ids[2*edge:2*edge+3]]
        line, _ = curve_samples(nodes, np.arange(3), tessellation_m)
        pieces.extend((line[:-1],))
        owner.extend([edge]*(len(line)-1))
    complete = np.vstack(pieces+[geometry[:, interface_ids[-1:]].T])
    start, edge_vector = complete[:-1], complete[1:]-complete[:-1]
    owner = np.asarray(owner)
    parameters = .5*(np.polynomial.legendre.leggauss(5)[0]+1.)
    samples = []; rejected_rays = 0
    def cross(a, b):
        return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]
    for edge in range((len(interface_ids)-1)//2):
        ids = interface_ids[2*edge:2*edge+3]
        vertices = geometry[:, ids]
        facet = facets[tuple(sorted((int(ids[0]), int(ids[-1]))))]
        adjacent = mesh.f2t[:, facet]
        cell = int(adjacent[adjacent >= 0][0])
        for q in parameters:
            position = vertices@np.array([(1-q)*(1-2*q), 4*q*(1-q), q*(2*q-1)])
            tangent = vertices@np.array([-3+4*q, 4-8*q, -1+4*q])
            # Purely diagnostic orientation classification, not a physical cutoff.
            if tangent[1] < abs(tangent[0]) or np.linalg.norm(tangent) <= 1e-30:
                continue
            outward = np.array([-tangent[1], tangent[0]])/np.linalg.norm(tangent)
            inward = -outward
            distance = (1.-position[0])/inward[0]
            wall_point = position+distance*inward
            if distance < 0. or wall_point[1] < -10. or wall_point[1] > SLOPE+1e-12:
                continue
            denominator = cross(inward, edge_vector)
            use = (abs(denominator) > 1e-24) & (owner != edge)
            ray_distance = np.full(len(start), np.inf)
            edge_parameter = np.zeros(len(start))
            ray_distance[use] = cross(start[use]-position, edge_vector[use])/denominator[use]
            edge_parameter[use] = cross(start[use]-position, inward)/denominator[use]
            crossings = use & (ray_distance > 2*tessellation_m) & (ray_distance < distance-2*tessellation_m) & (edge_parameter >= 0.) & (edge_parameter <= 1.)
            if np.any(crossings):
                rejected_rays += 1
                continue
            projection = outward@geometry[:, scalar.element_dofs[:, cell]]
            lower, upper = _quadratic_range(projection)
            width = upper-lower
            if width <= 0.:
                raise ValueError("nonpositive_element_normal_width")
            samples.append({"edge": edge, "cell": cell, "edge_parameter": float(q),
                "position_m": position.tolist(), "inward_normal": inward.tolist(),
                "right_wall_intersection_m": wall_point.tolist(), "gap_normal_m": float(distance),
                "adjacent_element_normal_width_m": float(width), "gap_over_element_normal_width": float(distance/width)})
    return {"scope": "DIAGNOSTIC_NOT_AUTOMATIC_THIN_REGION_ACCEPTANCE", "samples": len(samples),
        "wall_facing_orientation_rule": "upward tangent with |tau_z| >= |tau_x|; five interior Gauss points per full P2 edge",
        "rejected_rays_crossing_other_interface_before_wall": rejected_rays,
        "minimum_gap_witness": min(samples, key=lambda value: value["gap_normal_m"]) if samples else None,
        "minimum_gap_over_cell_width_witness": min(samples, key=lambda value: value["gap_over_element_normal_width"]) if samples else None,
        "crossing_test_tessellation_bound_m": tessellation_m,
        "element_normal_width": "exact quadratic projection range over the adjacent curved triangle",
        "interpretation": "A ratio near one indicates about one adjacent-cell normal width across the wall-connected resolved region. The gap necessarily tends to zero at P2; full-boundary h refinement, not this diagnostic alone, determines adequacy.",
        "molecular_retained_film_thickness_computed": False}
