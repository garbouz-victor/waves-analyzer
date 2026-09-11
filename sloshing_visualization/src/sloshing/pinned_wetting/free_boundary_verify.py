"""Independent PW2 material free-boundary audit.

This module does not import the production free-boundary solver, its residuals,
mesh motion, or diagnostics.  scikit-fem geometry/basis and SciPy assembly are a
shared numerical substrate, not an independently implemented FEM library.
"""
import hashlib
import json
import math
from pathlib import Path
from fractions import Fraction
from numpy.polynomial import Polynomial

import h5py
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, vstack
from skfem import (Basis, BilinearForm, ElementTriP1, ElementTriP2,
                   ElementVector, LinearForm, MeshTri2, asm)

MODEL_ID = "PW2_FULL_NS_LEFT_NAVIER_RIGHT_NOSLIP_MATERIAL_P2_V1"
METHOD = "full_material_P2_P1_nonlinear_midpoint"
SLOPE = float(np.tan(np.deg2rad(2.)))
E0 = 9.81*SLOPE*SLOPE/3.


def _hash_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def _canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def verifier_dependencies():
    """Cache identity includes independently implemented geometry helpers."""
    directory = Path(__file__).parent
    return {name: _hash_file(directory/name) for name in
            ("free_boundary_verify.py", "free_boundary_geometry_audit.py")}


def physical_contract():
    """Locked physics, independent of a possibly mutated native declaration."""
    return {
        "geometry": {"dimension": 2, "half_width_m": 1., "depth_m": 10.,
                     "mean_level_m": 0., "bottom_is_physical": True, "lid": False},
        "initial": {"interface": "z=x*tan(alpha)", "alpha_deg": 2.,
                    "velocity": "zero", "horizontal_forcing_after_release_m_s2": 0.},
        "material": {"kinematic_viscosity_m2_s": .01, "density_kg_m3": 1000.,
                     "gravity_m_s2": 9.81, "surface_tension_N_m": 0.,
                     "external_pressure": "constant_atmospheric", "gas_dynamics": "neglected"},
        "boundaries": {"left": {"type": "Navier", "impermeable": True, "slip_length_m": .5},
                       "right": {"type": "no_slip", "impermeable": True},
                       "bottom": {"type": "no_slip", "impermeable": True}},
        "wetting": {"irreversible": True, "connected_lower_closure": True,
                    "thickness_required": False, "evaporation": False, "dewetting": False,
                    "retained_layer": "inherited_zero_volume_one_way_state",
                    "resolved_film_must_be_mass_conserved": True, "initial_global_precursor": False},
        "right_contact": {"material_contact_at_P2_is_fixed": True,
                          "new_right_slip_forbidden": True,
                          "detached_bulk_height_is_not_wall_contact": True},
        "target_equations": "full_incompressible_Navier_Stokes_on_evolving_liquid_domain",
        "target_free_surface": "material_interface_with_actual_normals_and_zero_gauge_traction"}


@BilinearForm
def _scalar_mass(u, v, _):
    return u*v


@BilinearForm
def _strain(u, v, _):
    # nu*(grad u:grad v + grad u:(grad v)^T) = 2 nu D(u):D(v).
    return .01*sum(u.grad[i, j]*(v.grad[i, j]+v.grad[j, i])
                   for i in range(2) for j in range(2))


@BilinearForm
def _divergence(u, q, _):
    return q*(u.grad[0, 0]+u.grad[1, 1])


@BilinearForm
def _grad_div(u, v, _):
    return (u.grad[0, 0]+u.grad[1, 1])*(v.grad[0, 0]+v.grad[1, 1])


@LinearForm
def _grad_div_action(v, w):
    return w.gamma*w.dilation*(v.grad[0, 0]+v.grad[1, 1])


@LinearForm
def _one(v, _):
    return v


@BilinearForm
def _divergence_mass(u, v, w):
    return .5*w.dilation*u*v


@LinearForm
def _transfer_component(v, w):
    return w.previous_component*v


@LinearForm
def _momentum_action(v, w):
    """Action on stored fields; no production operator or sparse product."""
    acceleration = sum(w.acceleration[i]*v[i] for i in range(2))
    strain = .01*sum(w.state.grad[i, j]*(v.grad[i, j]+v.grad[j, i])
                     for i in range(2) for j in range(2))
    pressure = -w.pressure*(v.grad[0, 0]+v.grad[1, 1])
    geometric = .5*w.geometric_dilation*sum(w.state[i]*v[i] for i in range(2))
    grad_div = w.grad_div_gamma*(w.state.grad[0, 0]+w.state.grad[1, 1])*(v.grad[0, 0]+v.grad[1, 1])
    return acceleration+strain+pressure+geometric+grad_div


@LinearForm
def _weak_divergence_action(q, w):
    return q*w.dilation


def _scatter_scalar(matrix, components, n, component=None):
    sparse = matrix.tocoo()
    which = range(len(components)) if component is None else [component]
    row, col, data = [], [], []
    for k in which:
        row.extend(components[k][sparse.row]); col.extend(components[k][sparse.col])
        data.extend(sparse.data)
    return coo_matrix((data, (row, col)), shape=(n, n)).tocsr()


def _boundary_facets(mesh):
    boundary = mesh.boundary_facets()
    mids = mesh.p[:, mesh.facets[:, boundary]].mean(axis=1)
    checks = {"left": abs(mids[0]+1.) < 1e-12,
              "right": abs(mids[0]-1.) < 1e-12,
              "bottom": abs(mids[1]+10.) < 1e-12,
              "surface": abs(mids[1]-SLOPE*mids[0]) < 1e-12}
    return {key: boundary[mask] for key, mask in checks.items()}


def _ordered_interface(mesh, boundaries):
    scalar = Basis(mesh, ElementTriP2())
    links = {}
    for facet in boundaries["surface"]:
        a, b = map(int, mesh.facets[:, facet])
        midpoint = int(scalar.dofs.facet_dofs[0, facet])
        links.setdefault(a, []).append((b, midpoint))
        links.setdefault(b, []).append((a, midpoint))
    endpoints = [key for key, values in links.items() if len(values) == 1]
    if len(endpoints) != 2 or any(len(values) > 2 for values in links.values()):
        raise ValueError("initial_interface_is_not_one_complete_open_chain")
    here = min(endpoints, key=lambda key: mesh.p[0, key])
    ordered, visited = [here], set()
    while True:
        visited.add(here)
        choices = [(node, mid) for node, mid in links[here] if node not in visited]
        if not choices:
            break
        if len(choices) != 1:
            raise ValueError("ambiguous_interface_topology")
        here, mid = choices[0]
        ordered.extend((mid, here))
    if len(ordered) != 2*len(boundaries["surface"])+1:
        raise ValueError("missing_interface_branch")
    return np.asarray(ordered, dtype=int)


def _facet_quadrature(mesh, scalar, facets, intorder):
    """Reference-edge integration: no nonlinear global-to-cell inverse search.

    Curved-boundary FacetBasis may fail its generic invF Newton even when a
    valid material P2 cell is known. The facet's reference edge is exact native
    topology, so inversion is unnecessary. Normals use F^{-T} n_reference.
    """
    abscissa, weights = np.polynomial.legendre.leggauss(max(2, (intorder+2)//2))
    parameter, weights = .5*(abscissa+1.), .5*weights
    reference_vertices = np.array([[0., 0.], [1., 0.], [0., 1.]])
    barycentric_gradient = np.array([[-1., -1.], [1., 0.], [0., 1.]])
    for facet in facets:
        adjacent = mesh.f2t[:, facet]
        cell = int(adjacent[adjacent >= 0][0])
        endpoints = mesh.facets[:, facet]
        local = [int(np.flatnonzero(mesh.t[:, cell] == vertex)[0]) for vertex in endpoints]
        first, last = reference_vertices[local]
        edge_direction = last-first
        reference = first[:, None]+edge_direction[:, None]*parameter
        r, s = reference
        barycentric = np.array([1.-r-s, r, s])
        phi = np.empty((6, len(parameter)))
        gradient = np.empty((6, 2, len(parameter)))
        for i in range(3):
            phi[i] = barycentric[i]*(2*barycentric[i]-1.)
            gradient[i] = barycentric_gradient[i, :, None]*(4*barycentric[i]-1.)
        for i, (a, b) in enumerate(((0, 1), (1, 2), (0, 2)), start=3):
            phi[i] = 4*barycentric[a]*barycentric[b]
            gradient[i] = 4*(barycentric_gradient[a, :, None]*barycentric[b]
                             +barycentric_gradient[b, :, None]*barycentric[a])
        dofs = scalar.element_dofs[:, cell]
        nodes = mesh.p[:, dofs]
        coordinates = nodes@phi
        jacobian = np.einsum("ij,jkq->ikq", nodes, gradient)
        determinant = jacobian[0, 0]*jacobian[1, 1]-jacobian[0, 1]*jacobian[1, 0]
        if np.any(determinant == 0.):
            raise ValueError("zero_boundary_cell_Jacobian")
        inverse = np.array([[jacobian[1, 1], -jacobian[0, 1]],
                            [-jacobian[1, 0], jacobian[0, 0]]])/determinant
        reference_normal = np.array([edge_direction[1], -edge_direction[0]])
        opposite = reference_vertices[list({0, 1, 2}-set(local))[0]]
        if reference_normal@(opposite-.5*(first+last)) > 0.:
            reference_normal *= -1.
        normal = np.einsum("jiq,j->iq", inverse, reference_normal)
        normal /= np.linalg.norm(normal, axis=0)
        tangent = np.einsum("ijq,j->iq", jacobian, edge_direction)
        measure = weights*np.linalg.norm(tangent, axis=0)
        yield {"cell": cell, "dofs": dofs, "phi": phi, "gradient_reference": gradient,
               "barycentric": barycentric, "inverse": inverse, "coordinates": coordinates,
               "normal": normal, "measure": measure}


def _left_scalar_matrix(mesh, scalar, facets, intorder):
    rows, columns, values = [], [], []
    for q in _facet_quadrature(mesh, scalar, facets, intorder):
        local = .02*np.einsum("iq,jq,q->ij", q["phi"], q["phi"], q["measure"])
        rows.extend(np.repeat(q["dofs"], 6)); columns.extend(np.tile(q["dofs"], 6))
        values.extend(local.ravel())
    return coo_matrix((values, (rows, columns)), shape=(scalar.N, scalar.N)).tocsr()


def _free_surface_force(mesh, scalar, components, facets, intorder, nvelocity):
    force = np.zeros(nvelocity)
    for q in _facet_quadrature(mesh, scalar, facets, intorder):
        for axis in range(2):
            local = -9.81*np.einsum("iq,q,q,q->i", q["phi"], q["coordinates"][1], q["normal"][axis], q["measure"])
            np.add.at(force, components[axis][q["dofs"]], local)
    return force


def _traction_norms(op, v, p, side, intorder, hydrostatic, grad_div_gamma=0.):
    """Separate physical and numerical augmented natural-traction diagnostics.

    The isotropic numerical stress gamma*div(v)*I is not physical viscosity.
    On LEFT only the tangential Navier residual is appropriate; on FREE we
    report the complete physical, numerical and augmented traction vectors.
    Norms are quadrature diagnostics, not pointwise boundary/error certificates.
    """
    nodes = np.array([v[component] for component in op["components"]])
    squared = {name: 0. for name in ("physical", "grad_div", "augmented")}
    sampled_max = 0.
    for q in _facet_quadrature(op["mesh"], op["scalar"], op["mesh"].boundaries[side], intorder):
        velocity = nodes[:, q["dofs"]]@q["phi"]
        reference_gradient = np.einsum("ij,jkq->ikq", nodes[:, q["dofs"]], q["gradient_reference"])
        gradient = np.einsum("ijq,jkq->ikq", reference_gradient, q["inverse"])
        pressure = p[op["pressure"].element_dofs[:, q["cell"]]]@q["barycentric"]
        if hydrostatic:
            pressure -= 9.81*q["coordinates"][1]
        traction = .01*np.einsum("ijq,jq->iq", gradient+gradient.swapaxes(0, 1), q["normal"])-pressure*q["normal"]
        numerical = grad_div_gamma*(gradient[0, 0]+gradient[1, 1])*q["normal"]
        sampled_max = max(sampled_max, float(np.max(np.linalg.norm(numerical, axis=0))))
        if side == "left":
            values = {"physical": (traction[1]+.02*velocity[1])**2,
                      "grad_div": numerical[1]**2,
                      "augmented": (traction[1]+numerical[1]+.02*velocity[1])**2}
        else:
            values = {"physical": np.sum(traction**2, axis=0),
                      "grad_div": np.sum(numerical**2, axis=0),
                      "augmented": np.sum((traction+numerical)**2, axis=0)}
        for name, value in values.items():
            squared[name] += float(np.sum(value*q["measure"]))
    return {**{name: float(np.sqrt(value)) for name, value in squared.items()},
            "grad_div_sampled_max": sampled_max}


def _traction_norm(op, v, p, side, intorder, hydrostatic):
    return _traction_norms(op, v, p, side, intorder, hydrostatic)["physical"]


def _basis(geometry, triangles, boundaries, intorder=8):
    """Current-geometry spaces and asymmetric constraints, without operators."""
    mesh = MeshTri2(doflocs=np.asarray(geometry), t=np.asarray(triangles),
                    _boundaries=boundaries)
    scalar = Basis(mesh, ElementTriP2(), intorder=intorder)
    velocity = Basis(mesh, ElementVector(ElementTriP2()), intorder=intorder)
    pressure = Basis(mesh, ElementTriP1(), intorder=intorder)
    components = velocity.split_indices()
    wall = {name: velocity.get_dofs(name).all() for name in ("left", "right", "bottom")}
    fixed = np.unique(np.r_[np.intersect1d(wall["left"], components[0]),
                             wall["right"], wall["bottom"]])
    free = np.setdiff1d(np.arange(velocity.N), fixed)
    return {"mesh": mesh, "scalar": scalar, "velocity": velocity, "pressure": pressure,
            "components": components, "wall": wall, "fixed": fixed, "free": free}


def _loads(op, boundaries, intorder, hydrostatic_split):
    f_body = np.zeros(op["velocity"].N)
    f_body[op["components"][1]] = -9.81*asm(_one, op["scalar"])
    f = (_free_surface_force(op["mesh"], op["scalar"], op["components"],
                             boundaries["surface"], intorder, op["velocity"].N)
         if hydrostatic_split else f_body)
    return f, f_body


def rebuild(geometry, triangles, boundaries, intorder=8, *, hydrostatic_split=False,
            skew_divergence=False, advecting_velocity=None, grad_div_gamma=0.):
    """Full independent operators retained for startup and remap audits."""
    op = _basis(geometry, triangles, boundaries, intorder)
    mesh, scalar, velocity, pressure = (op[name] for name in ("mesh", "scalar", "velocity", "pressure"))
    components = op["components"]
    scalar_mass = asm(_scalar_mass, scalar)
    M = _scatter_scalar(scalar_mass, components, velocity.N)
    left_scalar = _left_scalar_matrix(mesh, scalar, boundaries["left"], intorder)
    left = _scatter_scalar(left_scalar, components, velocity.N, 1)
    f, f_body = _loads(op, boundaries, intorder, hydrostatic_split)
    geometric = csr_matrix(M.shape)
    if skew_divergence and advecting_velocity is not None:
        field = velocity.interpolate(advecting_velocity)
        dilation = field.grad[0, 0]+field.grad[1, 1]
        geometric = _scatter_scalar(asm(_divergence_mass, scalar, dilation=dilation), components, velocity.N)
    return {**op, "M": M, "bulk": asm(_strain, velocity).tocsr(), "left": left,
            "B": asm(_divergence, velocity, pressure).tocsr(), "f": f,
            "f_body": f_body, "geometric": geometric,
            "grad_div": (grad_div_gamma*asm(_grad_div, velocity)).tocsr() if grad_div_gamma else csr_matrix(M.shape)}


def _weak_divergence(op, field):
    return asm(_weak_divergence_action, op["pressure"],
               dilation=field.grad[0, 0]+field.grad[1, 1])


def _actions(op, velocity, acceleration, pressure, intorder, *,
             hydrostatic_split=False, skew_divergence=False, grad_div_gamma=0.):
    """Independent all-row midpoint actions and direct quadrature powers.

    Only the evaluation strategy differs from full sparse reconstruction.
    Initial unreduced Robin witnesses, startup and each mixed remap retain
    explicit matrices; no constrained row or physical term is inferred from
    production diagnostics.  Dissipation is evaluated directly from fields.
    """
    field = op["velocity"].interpolate(velocity)
    dilation = field.grad[0, 0]+field.grad[1, 1]
    residual = asm(_momentum_action, op["velocity"], state=field,
                   acceleration=op["velocity"].interpolate(acceleration),
                   pressure=op["pressure"].interpolate(pressure),
                   grad_div_gamma=grad_div_gamma,
                   geometric_dilation=dilation if skew_divergence else np.zeros_like(dilation))
    left_power = 0.
    for q in _facet_quadrature(op["mesh"], op["scalar"], op["mesh"].boundaries["left"], intorder):
        indices = op["components"][1][q["dofs"]]
        trace = velocity[indices]@q["phi"]
        local = .02*np.einsum("iq,q,q->i", q["phi"], trace, q["measure"])
        np.add.at(residual, indices, local)
        left_power += float(.02*np.sum(trace*trace*q["measure"]))
    f, f_body = _loads(op, op["mesh"].boundaries, intorder, hydrostatic_split)
    symmetric_gradient = field.grad+field.grad.swapaxes(0, 1)
    bulk_power = float(.005*np.sum(np.sum(symmetric_gradient**2, axis=(0, 1))*op["scalar"].dx))
    geometric_power = (float(.5*np.sum(dilation*np.sum(np.asarray(field)**2, axis=0)*op["scalar"].dx))
                       if skew_divergence else 0.)
    return {"momentum": residual-f, "weak_divergence": _weak_divergence(op, field),
            "field": field, "bulk_power": bulk_power, "left_power": left_power,
            "grad_div_power": float(grad_div_gamma*np.sum(dilation*dilation*op["scalar"].dx)),
            "geometric_power": geometric_power,
            "hydrostatic_power": float(velocity@(f-f_body))}


def determinant_minima(mesh, orientation):
    """Exact signed quadratic minima, using polynomial coefficient fitting.

    Independent implementation: fit [1,x,y,x²,xy,y²], enumerate all stationary
    points of that polynomial on a closed reference triangle. Initial cell
    orientation is supplied, never inferred from the potentially inverted cell.
    """
    pts = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
    monomials = np.array([[1, x, y, x*x, x*y, y*y] for x, y in pts])
    values = _determinant(mesh, pts.T)*np.asarray(orientation)[:, None]
    coefficients = values@np.linalg.inv(monomials).T
    minima = np.min(values[:, :3], axis=1)
    for cell, c in enumerate(coefficients):
        constant, x, y, xx, xy, yy = c
        candidates = []
        for start, direction in ((np.array([0., 0.]), np.array([1., 0.])),
                                 (np.array([0., 0.]), np.array([0., 1.])),
                                 (np.array([1., 0.]), np.array([-1., 1.]))):
            hessian = np.array([[2*xx, xy], [xy, 2*yy]])
            a = .5*direction@hessian@direction
            b = direction@(hessian@start+np.array([x, y]))
            if a != 0.:
                t = -b/(2*a)
                if 0. < t < 1.:
                    candidates.append(start+t*direction)
        hessian = np.array([[2*xx, xy], [xy, 2*yy]])
        if abs(np.linalg.det(hessian)) > 1e-28*max(np.linalg.norm(hessian)**2, 1e-30):
            stationary = np.linalg.solve(hessian, -np.array([x, y]))
            if np.min(stationary) > 0. and np.sum(stationary) < 1.:
                candidates.append(stationary)
        for a, b in candidates:
            minima[cell] = min(minima[cell], constant+x*a+y*b+xx*a*a+xy*a*b+yy*b*b)
    return minima


def _determinant(mesh, points):
    """No absolute value and no library exception before inversion evidence."""
    mapping = mesh.mapping()
    return mapping.J(0, 0, points)*mapping.J(1, 1, points)-mapping.J(0, 1, points)*mapping.J(1, 0, points)


def curve_samples(geometry, ids, tolerance_m=1e-6):
    """Whole ordered quadratic curve; analytic chord-error bound per segment."""
    curve, maximum_bound = [], 0.
    for k in range(0, len(ids)-2, 2):
        p = geometry[:, ids[k:k+3]]
        deviation = float(np.linalg.norm(p[:, 1]-.5*(p[:, 0]+p[:, 2])))
        count = max(1, int(np.ceil(np.sqrt(deviation/tolerance_m))))
        if count > 32768:
            raise ValueError("interface_tessellation_resource_bound")
        s = np.arange(count)/count
        basis = np.array([(1-s)*(1-2*s), 4*s*(1-s), s*(2*s-1)])
        curve.append((p@basis).T)
        maximum_bound = max(maximum_bound, deviation/count**2)
    curve.append(geometry[:, ids[-1:]].T)
    return np.vstack(curve), maximum_bound


def _curve_bounds(geometry, ids):
    minima = np.array([np.inf, np.inf]); maxima = -minima
    for k in range(0, len(ids)-2, 2):
        p = geometry[:, ids[k:k+3]]
        quadratic = 2*(p[:, 0]-2*p[:, 1]+p[:, 2])
        linear = -3*p[:, 0]+4*p[:, 1]-p[:, 2]
        for axis in range(2):
            values = [p[axis, 0], p[axis, 2]]
            if quadratic[axis] != 0.:
                q = -linear[axis]/(2*quadratic[axis])
                if 0. < q < 1.:
                    values.append(p[axis, 0]+q*linear[axis]+q*q*quadratic[axis])
            minima[axis] = min(minima[axis], *values)
            maxima[axis] = max(maxima[axis], *values)
    return minima, maxima


def _solid_contact_events(geometry, ids):
    """Exact quadratic wall-contact candidates, not a physical gap cutoff.

    Signed distances are positive inside the container.  The endpoint on the
    first LEFT edge and the endpoint on the last RIGHT edge are the only known
    contact roots.  Exact rational arithmetic represents the native IEEE
    coordinates and discriminates a double root from cancellation in float
    polynomial evaluation.  Strictly positive distances on the scale of eight
    input-coordinate ULPs are separately labelled precision-unresolved, not
    asserted to be physical intersections or a molecular-film thickness.
    """
    events = []
    for offset in range(0, len(ids)-2, 2):
        points = geometry[:, ids[offset:offset+3]]
        for wall, axis, level, sign in (("left", 0, -1., 1), ("right", 0, 1., -1),
                                         ("bottom", 1, -10., 1)):
            values = [sign*(Fraction.from_float(float(value))-Fraction.from_float(level))
                      for value in points[axis]]
            c, middle, end = values
            b = 4*middle-3*c-end
            a = 2*(c-2*middle+end)
            coincident = a == b == c == 0
            candidates = [(Fraction(1, 2), Fraction(0))] if coincident else [
                (Fraction(0), c), (Fraction(1), end)]
            if a > 0:
                stationary = -b/(2*a)
                if 0 < stationary < 1:
                    candidates.append((stationary, c-b*b/(4*a)))
            uncertainty = 8*max(math.ulp(float(value)) for value in points[axis])
            for parameter, distance in candidates:
                known = ((wall == "left" and offset == 0 and parameter == 0)
                         or (wall == "right" and offset == len(ids)-3 and parameter == 1))
                if known and distance == 0:
                    continue
                if distance > Fraction.from_float(uncertainty):
                    continue
                s = float(parameter)
                phi = np.array([(1-s)*(1-2*s), 4*s*(1-s), s*(2*s-1)])
                position = points@phi
                kind = ("coincident_surface_wall_segment" if coincident else
                        "crossing_or_outside" if distance < 0 else
                        "additional_exact_contact" if distance == 0 else
                        "precision_unresolved_contact")
                events.append({"wall": wall, "edge": offset//2, "edge_parameter": s,
                               "position_m": position.tolist(), "kind": kind,
                               "signed_distance_m": float(distance),
                               "input_coordinate_8ULP_scale_m": uncertainty})
    return events


def _left_height_is_global_record(heights, source):
    """Include every accepted endpoint, also an earlier unmarked plateau."""
    # max selects an existing native value: no interpolation, cancellation or
    # arithmetic tolerance is needed to decide the discrete history record.
    return float(heights[source]) == float(np.max(heights[:source+1]))


def _expected_rezone(controls, previous_time, last_rezone_time,
                     previous_exact_minima, initial_exact_minima):
    """Independent schedule from native geometry, never saved quality metadata.

    Inputs are per-cell signed minima already independently reconstructed on
    the entire closed quadratic triangles in the accepted-state geometry loop.
    Their ratio uses the fixed initial reference, not the most recent rezone.
    This dimensionless numerical adaptation is unrelated to physical slip.
    """
    interval = controls.get("rezone_interval_s", 0.)
    if interval <= 0.:
        return False
    scheduled = previous_time >= last_rezone_time+interval-1e-11
    adaptive = (controls.get("rezone_strategy", "harmonic_p2") == "quality_optimized"
                and float(np.min(np.asarray(previous_exact_minima)/initial_exact_minima)) < .01)
    return bool(scheduled or adaptive)


def _self_intersects(curve):
    def cross(a, b):
        return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]
    for i in range(len(curve)-3):
        a, r = curve[i], curve[i+1]-curve[i]
        b = curve[i+2:-1]; s = curve[i+3:]-b
        denominator = cross(r, s)
        use = abs(denominator) > 1e-24
        t = np.zeros(len(b)); u = t.copy()
        t[use] = cross(b[use]-a, s[use])/denominator[use]
        u[use] = cross(b[use]-a, r)/denominator[use]
        if np.any(use & (t > 1e-10) & (t < 1.-1e-10) & (u > 1e-10) & (u < 1.-1e-10)):
            return True
    return False


def _quadratic_curve_intersections(geometry, ids):
    """Residual-verified roots of complete P2 edges, not only their chords.

    Exact rational coefficients represent the native IEEE nodal coordinates.
    Eliminate one parameter to a degree-at-most-four polynomial, remove ONLY
    exact known shared-vertex factors, then find/polish numerical roots of the
    original two coordinate equations.  Collinear/coincident ranges are dealt
    with explicitly: they do not have isolated resultant roots.

    Positive transverse witnesses are concrete intersection evidence.  Finite
    precision roots, especially near tangencies, are not a rigorous separation
    certificate when the returned list is empty.  No physical gap or film-width
    cutoff is introduced.  Candidate/error tolerances concern algebra only.
    """
    geometry = np.asarray(geometry, dtype=float)
    if not np.isfinite(geometry[:, ids]).all():
        raise ValueError("nonfinite_quadratic_interface")
    def cross(a, b):
        return a[0]*b[1]-a[1]*b[0]
    def trim(p):
        p = list(p)
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p
    def add(p, q):
        return trim([(p[k] if k < len(p) else 0)+(q[k] if k < len(q) else 0)
                     for k in range(max(len(p), len(q)))])
    def times(p, q):
        result = [Fraction(0)]*(len(p)+len(q)-1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                result[i+j] += a*b
        return trim(result)
    def scaled(p, a):
        return trim([a*c for c in p])
    def at(p, x):
        value = Fraction(0)
        for c in p[::-1]:
            value = value*x+c
        return value
    def remove_exact_root(p, x):
        p = trim(p)
        while len(p) > 1 and at(p, x) == 0:
            quotient = [Fraction(0)]*(len(p)-1)
            quotient[-1] = p[-1]
            for k in range(len(p)-2, 0, -1):
                quotient[k-1] = p[k]+x*quotient[k]
            p = trim(quotient)
        return p
    def roots(p, shared=None):
        p = trim(p)
        if shared is not None:
            p = remove_exact_root(p, shared)
        result = []
        for endpoint in (Fraction(0), Fraction(1)):
            if len(p) > 1 and at(p, endpoint) == 0:
                result.append(float(endpoint))
                p = remove_exact_root(p, endpoint)
        if len(p) == 2:
            return result+[float(-p[0]/p[1])]
        if len(p) == 3:
            discriminant = p[1]*p[1]-4*p[0]*p[2]
            if discriminant < 0:
                return result  # exactly positive separation, not a near-real root
            if discriminant == 0:
                return result+[float(-p[1]/(2*p[2]))]
        if len(p) > 1:
            scale = max(abs(c) for c in p)
            for root in Polynomial([float(c/scale) for c in p]).roots():
                # Complex-near-real roots are seeds only, never evidence on
                # their own: both original coordinate equations follow below.
                if abs(root.imag) < 1e-6 and -.01 < root.real < 1.01:
                    result.append(float(root.real))
        return result
    def scalar_range(coefficients):
        c, b, a = coefficients
        values = [c, c+b+a]
        if a != 0 and 0 < -b/(2*a) < 1:
            values.append(c-b*b/(4*a))
        return min(values), max(values)

    edges = []; bounds = []
    for k in range(0, len(ids)-2, 2):
        p = geometry[:, ids[k:k+3]]
        r = [[Fraction.from_float(float(value)) for value in row] for row in p]
        exact = [[row[0], 4*row[1]-3*row[0]-row[2],
                  2*(row[0]-2*row[1]+row[2])] for row in r]
        edges.append(exact)
        control = np.array([p[:, 0], 2*p[:, 1]-.5*(p[:, 0]+p[:, 2]), p[:, 2]])
        # Outward padding protects broad-phase rejection from roundoff; it
        # neither changes the P2 curve nor decides whether there is contact.
        pad = 16*np.finfo(float).eps*max(float(np.max(abs(p))), np.finfo(float).tiny)
        bounds.append((control.min(axis=0)-pad, control.max(axis=0)+pad))
    events = []
    for first in range(len(edges)):
        for second in range(first+1, len(edges)):
            lo, hi = bounds[first]; other_lo, other_hi = bounds[second]
            if np.any(hi < other_lo) or np.any(other_hi < lo):
                continue
            adjacent = second == first+1
            pair_events = []
            for reverse in (False, True):
                i, j = (second, first) if reverse else (first, second)
                p, q = edges[i], edges[j]
                a, b, c = ([p[k][2] for k in (0, 1)], [p[k][1] for k in (0, 1)],
                           [p[k][0] for k in (0, 1)])
                d, e, f = ([q[k][2] for k in (0, 1)], [q[k][1] for k in (0, 1)],
                           [q[k][0] for k in (0, 1)])
                R = [[c[k]-f[k], b[k], a[k]] for k in (0, 1)]
                denominator = cross(e, d)
                shared = Fraction(0 if reverse else 1) if adjacent else None
                seeds = []; coincident = None
                if denominator != 0:
                    N = add(scaled(R[0], d[1]), scaled(R[1], -d[0]))
                    axis = int(abs(d[1]) > abs(d[0]))
                    equation = add(add(scaled(R[axis], denominator**2),
                                       scaled(N, -denominator*e[axis])),
                                   scaled(times(N, N), -d[axis]))
                    if all(coefficient == 0 for coefficient in equation):
                        tcoeff = [at(N, Fraction(0))/denominator,
                                  (N[1] if len(N) > 1 else 0)/denominator,
                                  (N[2] if len(N) > 2 else 0)/denominator]
                        lower, upper = scalar_range(tcoeff)
                        lower, upper = max(lower, Fraction(0)), min(upper, Fraction(1))
                        coincident = ("coincident_quadratic_overlap", lower, upper,
                                      Fraction(1 if reverse else 0) if adjacent else None)
                    else:
                        tpoly = Polynomial([float(coefficient/denominator) for coefficient in N])
                        seeds = [(s, float(tpoly(s))) for s in roots(equation, shared)]
                else:
                    direction = e if max(map(abs, e)) >= max(map(abs, d)) else d
                    if direction == [0, 0]:
                        pair_events.append({"kind": "collapsed_quadratic_edge", "edge": j})
                        continue
                    axis = int(abs(direction[1]) > abs(direction[0]))
                    equation = add(scaled(R[0], direction[1]), scaled(R[1], -direction[0]))
                    if all(coefficient == 0 for coefficient in equation):
                        l0, h0 = scalar_range(p[axis]); l1, h1 = scalar_range(q[axis])
                        coincident = ("collinear_quadratic_overlap", max(l0, l1), min(h0, h1),
                                      Fraction.from_float(float(geometry[axis, ids[2*second]])) if adjacent else None)
                    else:
                        for s in roots(equation, shared):
                            value = float(R[axis][0])+float(R[axis][1])*s+float(R[axis][2])*s*s
                            for t in Polynomial([-value, float(e[axis]), float(d[axis])]).roots():
                                if abs(t.imag) < 1e-6:
                                    seeds.append((s, float(t.real)))
                if coincident is not None:
                    kind, lower, upper, known = coincident
                    if lower < upper or (lower == upper and lower != known):
                        pair_events.append({"kind": kind if lower < upper else "nonincident_exact_contact",
                                            "exact_overlap_interval": [str(lower), str(upper)]})
                    continue
                A, B, D, E = map(lambda values: np.asarray(values, dtype=float), (a, b, d, e))
                offset = np.asarray([float(c[k]-f[k]) for k in (0, 1)])
                coordinate_scale = max(float(np.max(abs(geometry[:, ids[2*first:2*second+3]]))), 1e-30)
                residual_tolerance = 64*np.finfo(float).eps*coordinate_scale
                for seed in seeds:
                    parameter = np.asarray(seed, dtype=float)
                    if not np.isfinite(parameter).all() or np.min(parameter) < -.01 or np.max(parameter) > 1.01:
                        continue
                    for _ in range(25):
                        s, t = parameter
                        residual = offset+B*s+A*s*s-E*t-D*t*t
                        derivative = np.column_stack((B+2*A*s, -E-2*D*t))
                        try:
                            change = np.linalg.solve(derivative, residual)
                        except np.linalg.LinAlgError:
                            break
                        parameter -= change
                        if np.linalg.norm(change, np.inf) <= 8*np.finfo(float).eps:
                            break
                    s, t = parameter
                    residual = float(np.linalg.norm(offset+B*s+A*s*s-E*t-D*t*t))
                    if not (np.isfinite(parameter).all() and residual <= residual_tolerance
                            and np.min(parameter) >= 0. and np.max(parameter) <= 1.):
                        continue
                    ordered = parameter[::-1] if reverse else parameter
                    # Known factors have already been removed exactly.  This
                    # last equality only excludes a numerical iteration that
                    # landed on the actual shared root, not a nearby interval.
                    if adjacent and np.array_equal(ordered, [1., 0.]):
                        continue
                    if any("parameters" in entry and np.max(abs(np.asarray(entry["parameters"])-ordered))
                           <= 64*np.finfo(float).eps for entry in pair_events):
                        continue
                    tangent0, tangent1 = B+2*A*s, E+2*D*t
                    norm_product = float(np.linalg.norm(tangent0)*np.linalg.norm(tangent1))
                    normalized_cross = float(cross(tangent0, tangent1)/norm_product) if norm_product else 0.
                    transverse = bool(abs(normalized_cross) > 1024*np.finfo(float).eps)
                    exact_contact = all(at(p[k], Fraction.from_float(float(s))) ==
                                        at(q[k], Fraction.from_float(float(t))) for k in (0, 1))
                    kind = ("residual_verified_quadratic_intersection" if transverse else
                            "exact_quadratic_contact" if exact_contact else
                            "precision_unresolved_quadratic_contact")
                    pair_events.append({"kind": kind,
                        "parameters": ordered.tolist(),
                        "position_m": (np.asarray(c, dtype=float)+B*s+A*s*s).tolist(),
                        "equation_residual_m": residual,
                        "algebraic_residual_tolerance_m": residual_tolerance,
                        "strictly_interior": bool(np.min(ordered) > 0. and np.max(ordered) < 1.),
                        "normalized_tangent_cross": normalized_cross,
                        "transverse": transverse,
                        "exact_rational_coordinate_equality_at_float_parameters": exact_contact})
            # The reverse elimination is an additional numerical conditioning
            # check, not a distinct geometric event for coincident ranges.
            for event in pair_events:
                witness = {"edges": [first, second], **event}
                if witness not in events:
                    events.append(witness)
    return events


def _measure(op, v):
    field = np.asarray(op["velocity"].interpolate(v))
    dx = op["scalar"].dx
    z = np.asarray(op["scalar"].global_coordinates())[1]
    kinetic = .5*float(np.sum(np.sum(field*field, axis=0)*dx))
    potential = 9.81*(float(np.sum(z*dx))+100.)
    return {"volume_m2": float(np.sum(dx)), "kinetic": kinetic,
            "potential": potential, "energy": kinetic+potential}


def _read_matrix(group):
    return csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]),
                      shape=tuple(group.attrs["shape"]))


def _audit_rezone(group, old, old_v, triangles, boundaries, controls, ids, orientation):
    """Read-only transfer audit, no production search/projection implementation."""
    new_X, new_v = group["geometry"][:], group["velocity"][:]
    new = rebuild(new_X, triangles, boundaries, controls["intorder"],
                  hydrostatic_split=controls.get("hydrostatic_split", False))
    metrics = {}; errors = []
    def require(ok, reason):
        if not ok:
            errors.append(reason)
    minimum = determinant_minima(new["mesh"], orientation)
    require(np.min(minimum) > 0., "rezone_invalid_curved_geometry")
    metrics["rezone_surface_geometry_change_m"] = float(np.max(abs(new_X[:, ids]-old["mesh"].p[:, ids])))
    metrics["rezone_constrained_velocity_m_s"] = float(np.max(abs(new_v[new["fixed"]])))
    for side in ("left", "right", "bottom"):
        nodes = new["scalar"].get_dofs(side).all()
        axis, value = (1, -10.) if side == "bottom" else (0, -1. if side == "left" else 1.)
        require(np.max(abs(new_X[axis, nodes]-value)) <= 1e-12, "rezone_changed_solid_wall_geometry")
    cells, reference = group["source_cells"][:], group["source_reference_coordinates"][:]
    query = np.asarray(new["scalar"].global_coordinates())
    require(cells.shape == (query.shape[1]*query.shape[2],) and reference.shape == (2, len(cells)),
            "rezone_inverse_map_shape")
    require(np.issubdtype(cells.dtype, np.integer) and np.all(cells >= 0) and np.all(cells < triangles.shape[1]),
            "rezone_wrong_source_cells")
    require(np.isfinite(reference).all() and np.min(reference) >= -1e-10
            and np.max(reference.sum(axis=0)) <= 1.+1e-10, "rezone_extrapolated_outside_old_liquid")
    if errors:
        raise ValueError(",".join(errors))
    # Barycentric nodal interpolation, independent of production inverse/search.
    r, s = reference
    l0, l1, l2 = 1.-r-s, r, s
    basis = np.array([l0*(2*l0-1), l1*(2*l1-1), l2*(2*l2-1),
                      4*l0*l1, 4*l1*l2, 4*l0*l2])
    dofs = old["scalar"].element_dofs[:, cells]
    mapped = np.einsum("ijn,jn->in", old["mesh"].p[:, dofs], basis)
    metrics["rezone_inverse_map_error_m"] = float(np.max(np.linalg.norm(mapped-query.reshape(2, -1), axis=0)))
    old_nodes = np.array([old_v[component] for component in old["components"]])
    values = np.einsum("ijn,jn->in", old_nodes[:, dofs], basis).reshape(query.shape)
    rhs = np.zeros(new["velocity"].N)
    P = np.zeros((2, new["velocity"].N))
    for axis, component in enumerate(new["components"]):
        rhs[component] = asm(_transfer_component, new["scalar"], previous_component=values[axis])
        P[axis, component] = asm(_one, new["scalar"])
    old_field = np.asarray(old["velocity"].interpolate(old_v))
    target = np.sum(old_field*old["scalar"].dx[None, :, :], axis=(1, 2))
    free = new["free"]
    C = vstack((new["B"][:, free], csr_matrix(P[:, free])), format="csr")
    multipliers = group["projection_multipliers"][:]
    equation = new["M"][free][:, free]@new_v[free]+C.T@multipliers-rhs[free]
    constraint = C@new_v[free]-np.r_[np.zeros(new["B"].shape[0]), target]
    metrics["rezone_projection_residual"] = float(max(np.max(abs(equation)), np.max(abs(constraint))))
    metrics["rezone_momentum_error"] = float(np.max(abs(P@new_v-target)))
    metrics["rezone_weak_divergence"] = float(np.max(abs(new["B"]@new_v)))
    new_field = np.asarray(new["velocity"].interpolate(new_v))
    norm_error = float(np.sum(np.sum((new_field-values)**2, axis=0)*new["scalar"].dx))
    norm_old = float(np.sum(np.sum(values**2, axis=0)*new["scalar"].dx))
    metrics["rezone_L2_velocity_projection_relative"] = float(np.sqrt(norm_error/max(norm_old, 1e-30)))
    surface_velocity_change = np.array([new_v[component][ids]-old_v[old["components"][axis]][ids]
                                        for axis, component in enumerate(new["components"])])
    metrics["rezone_surface_velocity_change_m_s"] = float(np.max(abs(surface_velocity_change)))
    old_measure, new_measure = _measure(old, old_v), _measure(new, new_v)
    energy_work = new_measure["energy"]-old_measure["energy"]
    volume_change = new_measure["volume_m2"]-old_measure["volume_m2"]
    metrics["rezone_volume_change_m2"] = abs(volume_change)
    metrics["rezone_saved_energy_work_discrepancy"] = abs(float(group.attrs["energy_work"])-energy_work)
    metrics["rezone_saved_volume_change_discrepancy"] = abs(float(group.attrs["volume_change"])-volume_change)
    return new, new_v, energy_work, metrics, errors


def verify_native(native, *, require_target=False, heartbeat=lambda: None):
    """Audit every committed state; short-run PASS is explicitly not 5s release."""
    path = Path(native)
    errors = []
    maxima = {}
    witnesses = {}
    current_index = None
    last_complete = None
    report = {"passed": False, "scope": "INDEPENDENT_FULL_MATERIAL_NS_NATIVE",
              "model_id": MODEL_ID, "native": str(path), "require_target": require_target,
              "shared_substrate": "scikit-fem geometry/basis and scipy; no production operators/residuals",
              "operator_evaluation": "independent all-row LinearForm actions and quadrature powers; full sparse startup, unreduced initial LEFT-only witness and every remap",
              "verifier_sha256": _hash_file(__file__), "verifier_dependencies": verifier_dependencies(),
              "audit_complete": False, "additional_solid_contact_events": [],
              "quadratic_interface_intersection_events": [],
              "quadratic_curve_configurations_checked": 0,
              "quadratic_intersection_method": "exact native-coordinate polynomial coefficients and shared-root factor removal, numerical elimination roots and original-P2 residual verification; collinear/coincident ranges checked separately",
              "quadratic_intersection_qualification": "Every accepted and actual collocation curve is checked independently of production chords. Empty finite-precision root sets are not rigorous tangency/separation certificates. Residual-verified transverse roots establish concrete intersections; near-tangent algebraic witnesses retain numerical uncertainty. No physical gap cutoff is used."}
    def require(ok, name):
        if not ok and name not in errors:
            errors.append(name)
    def observe(name, value, step):
        value = float(value)
        if not np.isfinite(value):
            require(False, "nonfinite:"+name)
            return
        if name not in maxima or value > maxima[name]:
            maxima[name] = value
            witnesses[name] = int(step)
    def check_quadratic(geometry, ids, step, time_s, midpoint=False):
        events = _quadratic_curve_intersections(geometry, ids)
        report["quadratic_curve_configurations_checked"] += 1
        for event in events:
            report["quadratic_interface_intersection_events"].append({
                "step": int(step), "time_s": float(time_s),
                "configuration": "midpoint_collocation" if midpoint else "accepted_endpoint", **event})
        prefix = "midpoint_" if midpoint else ""
        unresolved = [event for event in events if event["kind"] == "precision_unresolved_quadratic_contact"]
        require(not unresolved, prefix+"precision_unresolved_full_quadratic_contact")
        require(len(unresolved) == len(events), prefix+"full_quadratic_interface_intersection_or_overlap")
    try:
        with h5py.File(path, "r") as h:
            ident = json.loads(h.attrs["numerical_identity"])
            require(h.attrs["model_id"] == MODEL_ID and ident.get("model_id") == MODEL_ID,
                    "not_new_full_nonlinear_model")
            require(h.attrs["equations"] == "full_material_NS" and ident.get("method") == METHOD,
                    "wrong_equations_or_ALE_method")
            require(not bool(h.attrs["synthetic"]) and ident.get("synthetic") is False,
                    "synthetic_or_non_solver_source")
            require(ident.get("physical_contract") == physical_contract(), "locked_physical_contract_mismatch")
            require(h.attrs["numerical_digest"] == _canonical(ident), "numerical_identity_digest_mismatch")
            require(json.loads((path.parent/"identity.json").read_text()) == ident,
                    "identity_sidecar_mismatch")
            resolved = json.loads((path.parent/"resolved_case.json").read_text())
            require(resolved.get("physical_contract") == physical_contract(), "resolved_case_physics_mismatch")
            for name, digest in ident["source_hashes"].items():
                require(_hash_file(path.parent/"source_snapshot"/name) == digest,
                        "numerical_source_snapshot_digest_mismatch:"+name)
            require("pinned_wetting/free_boundary.py" in ident["source_hashes"], "not_new_solver_snapshot")
            controls = ident["controls"]
            hydrostatic = controls.get("hydrostatic_split", False)
            skew_divergence = controls.get("skew_divergence", False)
            gamma = controls.get("grad_div_gamma_m2_s", 0.)
            require(isinstance(gamma, (int, float)) and np.isfinite(gamma) and gamma >= 0.,
                    "invalid_numerical_grad_div_coefficient")
            declared_grad_div = "grad_div_gamma_m2_s" in controls
            pressure_variable = "pi=p+gz" if hydrostatic else "gauge_p"
            require(h.attrs.get("pressure_variable", "gauge_p") == pressure_variable,
                    "pressure_variable_does_not_match_equations")
            topology = h["topology"]
            triangles, initial = topology["triangles"][:], topology["initial_geometry"][:]
            bare = MeshTri2(doflocs=initial, t=triangles)
            boundaries = _boundary_facets(bare)
            require(sum(map(len, boundaries.values())) == len(bare.boundary_facets()),
                    "initial_boundary_not_exact_container_and_line")
            for side, facets in boundaries.items():
                require(np.array_equal(np.sort(facets), np.sort(topology[side+"_facets"][:])),
                        "wrong_boundary_facets_or_swapped_sides:"+side)
            ids = _ordered_interface(bare, boundaries)
            require(np.array_equal(topology["interface_nodes"][:], ids), "missing_or_reordered_interface_branch")
            op0 = rebuild(initial, triangles, boundaries, controls["intorder"],
                          hydrostatic_split=hydrostatic, grad_div_gamma=gamma)
            orientations = np.sign(_determinant(bare, np.array([[0.], [0.]]))[:, 0])
            require(np.array_equal(topology["orientation"][:], orientations), "wrong_initial_orientation")
            initial_minimum = determinant_minima(bare, orientations)
            require(np.array_equal(topology["fixed_velocity_dofs"][:], op0["fixed"]),
                    "right_freed_or_wrong_velocity_constraints")
            require(np.array_equal(topology["component_dofs"][:], np.asarray(op0["components"])),
                    "wrong_normal_tangential_components")
            saved_left = _read_matrix(h["initial_left_wall_operator"])
            expected_left = op0["left"]
            difference = saved_left-expected_left
            scale = max(float(np.max(abs(expected_left.data))), 1e-30)
            require(not difference.nnz or np.max(abs(difference.data)) <= 1e-12*scale,
                    "full_LEFT_Robin_matrix_mismatch")
            right = op0["wall"]["right"]
            forbidden = saved_left[right]
            require(not forbidden.nnz or np.max(abs(forbidden.data)) <= 1e-14*scale,
                    "forbidden_RIGHT_Robin_even_on_constrained_DOFs")
            require(not declared_grad_div or "initial_grad_div_operator" in h,
                    "missing_initial_grad_div_operator")
            if "initial_grad_div_operator" in h:
                saved_g = _read_matrix(h["initial_grad_div_operator"])
                expected_g = op0["grad_div"]
                gscale = max(float(np.max(abs(expected_g.data))) if expected_g.nnz else 0., 1.)
                delta = saved_g-expected_g
                require(not delta.nnz or np.max(abs(delta.data)) <= 1e-12*gscale,
                        "initial_grad_div_operator_mismatch")
                symmetry = saved_g-saved_g.T
                require(not symmetry.nnz or np.max(abs(symmetry.data)) <= 1e-12*gscale,
                        "initial_grad_div_operator_not_symmetric")
                # Independent Gram-form equality establishes the PSD form.
                # These full-DOF witnesses additionally test action/power; they
                # are not mislabelled as an exhaustive numerical eigenspectrum.
                energies = []; action_errors = []
                indices = np.arange(op0["velocity"].N, dtype=float)
                for values in (np.sin(.31*indices), np.cos(.17*indices), np.ones_like(indices)):
                    field = op0["velocity"].interpolate(values)
                    dilation = field.grad[0, 0]+field.grad[1, 1]
                    independent_action = asm(_grad_div_action, op0["velocity"], gamma=gamma, dilation=dilation)
                    actual_action = saved_g@values
                    action_error = float(np.max(abs(actual_action-independent_action))
                                         /max(1., float(np.max(abs(independent_action)))))
                    energy = float(values@actual_action)
                    independent_power = float(gamma*np.sum(dilation*dilation*op0["scalar"].dx))
                    require(action_error <= 1e-11, "initial_grad_div_action_mismatch")
                    require(energy >= -1e-11*gscale and independent_power >= -1e-12,
                            "initial_grad_div_negative_quadratic_form")
                    require(abs(energy-independent_power) <= 1e-10*max(1., abs(independent_power)),
                            "initial_grad_div_quadratic_power_mismatch")
                    energies.append(energy); action_errors.append(action_error)
                report["initial_grad_div_witness"] = {"gamma_m2_s": gamma,
                    "full_matrix_matches_independent_divergence_Gram_form": not delta.nnz or bool(np.max(abs(delta.data)) <= 1e-12*gscale),
                    "quadratic_form_witnesses": energies, "action_relative_errors": action_errors,
                    "PSD_basis": "nonnegative gamma times integral of squared divergence; matrix equality and deterministic all-DOF action witnesses, not an eigenspectrum claim"}
            # Natural side normals are independently geometric, not metadata.
            for side, normal in (("left", [-1., 0.]), ("right", [1., 0.]), ("bottom", [0., -1.])):
                actual = op0["velocity"].boundary(side).normals
                require(np.max(abs(actual-np.asarray(normal)[:, None, None])) <= 1e-12,
                        "wrong_outward_normal:"+side)
            require(np.max(abs(initial[1, ids]-SLOPE*initial[0, ids])) <= 1e-12,
                    "initial_surface_not_prescribed_straight_line")
            require(np.max(abs(initial[:, ids[[0, -1]]]-np.array([[-1., 1.], [-SLOPE, SLOPE]]))) <= 1e-12,
                    "wrong_initial_P1_P2")
            keys = sorted(h["states"], key=int)
            require([int(key) for key in keys] == list(range(len(keys))), "missing_or_duplicate_accepted_step")
            require(int(h.attrs["last_accepted_step"]) == len(keys)-1, "uncommitted_states_in_accepted_history")
            run_id = h.attrs["run_id"]
            contact_history, time_history, previous_markers = [], [], []
            H = np.array([-SLOPE, SLOPE])
            cumulative_bulk = cumulative_left = cumulative_kinetic = cumulative_potential = 0.
            cumulative_rezone = 0.; last_rezone_time = 0.; rezone_count = 0
            cumulative_grad_div = 0.
            previous_X = previous_v = previous_op = previous_energy = None
            previous_minimum = None
            max_tessellation = 0.
            for index, key in enumerate(keys):
                current_index = index
                group = h["states"][key]
                X, v, p = (group[name][:] for name in ("geometry", "velocity", "pressure_midpoint"))
                require(bool(group.attrs["accepted"]) and int(group.attrs["step"]) == index,
                        "rejected_or_wrong_step_in_accepted_history")
                require(np.isfinite(X).all() and np.isfinite(v).all() and np.isfinite(p).all(), "nonfinite_fields")
                time_s = float(group.attrs["time_s"])
                require(time_s == 0. if index == 0 else time_s > time_history[-1], "bad_accepted_times")
                op = op0 if index == 0 else _basis(X, triangles, boundaries, controls["intorder"])
                minimum = determinant_minima(op["mesh"], orientations)
                require(np.min(minimum) > 0., "nonpositive_curved_cell_Jacobian")
                observe("inverse_min_Jacobian", 1./max(float(np.min(minimum)), 1e-300), index)
                scalar_values = np.asarray([v[component] for component in op["components"]])
                observe("right_velocity_max_m_s", np.max(abs(v[op["wall"]["right"]])), index)
                observe("bottom_velocity_max_m_s", np.max(abs(v[op["wall"]["bottom"]])), index)
                left_normal = np.intersect1d(op["wall"]["left"], op["components"][0])
                observe("left_normal_velocity_max_m_s", np.max(abs(v[left_normal])), index)
                observe("accepted_endpoint_weak_divergence", np.max(abs(_weak_divergence(op, op["velocity"].interpolate(v)))), index)
                observe("P2_position_drift_m", np.max(abs(X[:, ids[-1]]-initial[:, ids[-1]])), index)
                for side in ("left", "right", "bottom"):
                    nodes = op["scalar"].get_dofs(side).all()
                    axis, value = (1, -10.) if side == "bottom" else (0, -1. if side == "left" else 1.)
                    observe("solid_boundary_position_drift_m", np.max(abs(X[axis, nodes]-value)), index)
                bounds_min, bounds_max = _curve_bounds(X, ids)
                require(bounds_min[0] >= -1.-1e-12 and bounds_max[0] <= 1.+1e-12
                        and bounds_min[1] > -10., "full_interface_crosses_solid_wall")
                for event in _solid_contact_events(X, ids):
                    report["additional_solid_contact_events"].append({"step": index, "time_s": time_s, **event})
                    require(False, "precision_unresolved_extra_solid_contact" if
                            event["kind"] == "precision_unresolved_contact" else "unsupported_extra_solid_contact")
                curve, tess = curve_samples(X, ids)
                max_tessellation = max(max_tessellation, tess)
                require(not _self_intersects(curve), "interface_self_intersection_at_controlled_tessellation")
                check_quadratic(X, ids, index, time_s)
                measured = _measure(op, v)
                observe("volume_relative", abs(measured["volume_m2"]-20.)/20., index)
                row = json.loads(group.attrs["diagnostics"])
                if declared_grad_div:
                    require("numerical_grad_div_work" in row and "cumulative_numerical_grad_div_work" in row,
                            "missing_numerical_grad_div_work")
                require(row.get("numerical_grad_div_work", 0.) <= 0.
                        and row.get("cumulative_numerical_grad_div_work", 0.) <= 0.,
                        "numerical_grad_div_work_wrong_sign")
                for name in ("volume_m2", "kinetic", "potential", "energy"):
                    observe("saved_diagnostic_discrepancy", abs(row[name]-measured[name]), index)
                R = X[1, ids[[0, -1]]]
                H = np.maximum(H, R)
                coating = json.loads(group.attrs["coating"])
                require(set(coating) == {"run_id", "bottom", "H", "markers", "accepted_step", "time"},
                        "retained_state_not_single_lower_connected_intervals")
                require(coating["bottom"] == -10. and coating["run_id"] == run_id,
                        "retained_state_wrong_container_or_run")
                require(coating["accepted_step"] == index and coating["time"] == time_s,
                        "coating_history_lost_on_restart_or_rejected_trial")
                observe("retained_running_max_error_m", np.max(abs(np.asarray(coating["H"])-H)), index)
                observe("right_H_error_m", abs(coating["H"][1]-SLOPE), index)
                contact_history.append(R.tolist()); time_history.append(time_s)
                markers = coating["markers"]
                require(markers[:len(previous_markers)] == previous_markers, "historical_marker_mutation_or_history_loss")
                require(len({m["marker_id"] for m in markers}) == len(markers), "duplicate_marker_ID")
                initial_markers = [
                    {"marker_id": "P1", "side": "L", "height_m": float(initial[1, ids[0]]), "created_at_s": 0., "state_id": f"{run_id}:0"},
                    {"marker_id": "P2", "side": "R", "height_m": float(initial[1, ids[-1]]), "created_at_s": 0., "state_id": f"{run_id}:0"}]
                require(markers[:2] == initial_markers, "P1_P2_not_immutable_initial_coordinates")
                for marker in markers[2:]:
                    require(marker["side"] == "L", "new_RIGHT_record_forbidden")
                    source = int(marker["state_id"].rsplit(":", 1)[1])
                    require(marker["state_id"].startswith(str(run_id)+":"), "marker_wrong_native_run")
                    if 0 < source < index:
                        heights = np.asarray(contact_history)[:, 0]
                        require(abs(marker["height_m"]-heights[source]) <= 1e-12
                                and heights[source] > heights[source-1]
                                and heights[source] > heights[source+1]
                                and marker["created_at_s"] == time_history[source+1],
                                "marker_not_actual_confirmed_LEFT_peak")
                        require(_left_height_is_global_record(heights, source),
                                "marker_not_global_LEFT_record_at_source")
                    else:
                        require(False, "marker_not_confirmed_accepted_source_step")
                expected_new = []
                if index >= 2:
                    heights = np.asarray(contact_history)[:, 0]
                    prior_record = max(marker["height_m"] for marker in previous_markers if marker["side"] == "L")
                    if (heights[-3] < heights[-2] > heights[-1] and heights[-2] > prior_record+1e-6
                            and _left_height_is_global_record(heights, len(heights)-2)):
                        expected_new = [{"marker_id": f"P{len(previous_markers)+1}", "side": "L",
                                         "height_m": float(heights[-2]), "created_at_s": time_s,
                                         "state_id": f"{run_id}:{index-1}"}]
                if index:
                    require(markers[len(previous_markers):] == expected_new, "record_detector_missing_or_invented_LEFT_peak")
                previous_markers = markers
                if index == 0:
                    require(np.array_equal(X, initial) and np.max(abs(v)) == 0., "initial_geometry_or_velocity_changed")
                    acceleration = h["initial_acceleration"][:]
                    startup = op["M"]@acceleration-op["B"].T@p-op["f"]
                    observe("startup_momentum", np.max(abs(startup[op["free"]])), index)
                    observe("startup_weak_divergence", np.max(abs(op["B"]@acceleration)), index)
                    observe("startup_fixed_acceleration", np.max(abs(acceleration[op["fixed"]])), index)
                    observe("saved_grad_div_work_discrepancy", abs(row.get("numerical_grad_div_work", 0.)), index)
                else:
                    dt = time_s-time_history[-2]
                    require(abs(row["dt_s"]-dt) <= 1e-12, "saved_dt_not_physical_step_interval")
                    start_X, start_v, start_op = previous_X, previous_v, previous_op
                    start_energy = previous_energy
                    expected_rezone = _expected_rezone(controls, time_history[-2], last_rezone_time,
                                                       previous_minimum, initial_minimum)
                    require(("rezoning_before_step" in group) == expected_rezone,
                            "rezone_schedule_or_native_transfer_evidence_missing")
                    rezone_work = 0.
                    if "rezoning_before_step" in group:
                        start_op, start_v, rezone_work, transfer_metrics, transfer_errors = _audit_rezone(
                            group["rezoning_before_step"], previous_op, previous_v,
                            triangles, boundaries, controls, ids, orientations)
                        for name, value in transfer_metrics.items():
                            observe(name, value, index)
                        for name in transfer_errors:
                            require(False, name)
                        start_X = start_op["mesh"].p
                        start_energy = _measure(start_op, start_v)["energy"]
                        cumulative_rezone += rezone_work
                        last_rezone_time = time_history[-2]
                        rezone_count += 1
                    observe("saved_rezone_work_discrepancy", abs(row.get("rezone_work", 0.)-rezone_work), index)
                    require(abs(row.get("last_rezone_time_s", 0.)-last_rezone_time) <= 1e-12,
                            "rezone_history_lost_on_restart")
                    vm = .5*(v+start_v); dv = v-start_v
                    midpoint = _basis(.5*(X+start_X), triangles, boundaries, controls["intorder"])
                    nodal = np.asarray([vm[component] for component in midpoint["components"]])
                    observe("material_kinematic_defect_m", np.max(abs(X-start_X-dt*nodal)), index)
                    action = _actions(midpoint, vm, dv/dt, p, controls["intorder"],
                                      hydrostatic_split=hydrostatic, skew_divergence=skew_divergence,
                                      grad_div_gamma=gamma)
                    observe("momentum_residual", np.max(abs(action["momentum"][midpoint["free"]])), index)
                    observe("weak_divergence", np.max(abs(action["weak_divergence"])), index)
                    min_mid = determinant_minima(midpoint["mesh"], orientations)
                    require(np.min(min_mid) > 0., "midpoint_geometry_invalid")
                    midpoint_X = midpoint["mesh"].p
                    midpoint_curve, midpoint_tess = curve_samples(midpoint_X, ids)
                    max_tessellation = max(max_tessellation, midpoint_tess)
                    require(not _self_intersects(midpoint_curve),
                            "midpoint_interface_self_intersection_at_controlled_tessellation")
                    check_quadratic(midpoint_X, ids, index, .5*(time_s+time_history[-2]), midpoint=True)
                    for event in _solid_contact_events(midpoint_X, ids):
                        report["additional_solid_contact_events"].append({
                            "step": index, "time_s": .5*(time_s+time_history[-2]),
                            "configuration": "midpoint_collocation", **event})
                        require(False, "midpoint_precision_unresolved_extra_solid_contact" if
                                event["kind"] == "precision_unresolved_contact" else
                                "midpoint_unsupported_extra_solid_contact")
                    field = action["field"]
                    ref = midpoint["scalar"].X
                    J0 = _determinant(start_op["mesh"], ref)
                    J1 = _determinant(op["mesh"], ref)
                    Jm = _determinant(midpoint["mesh"], ref)
                    divergence = field.grad[0, 0]+field.grad[1, 1]
                    observe("stage_strong_divergence_Linf_s_inv", np.max(abs(divergence)), index)
                    observe("stage_strong_divergence_L2", np.sqrt(np.sum(divergence**2*midpoint["scalar"].dx)), index)
                    observe("maximum_step_dilation", dt*np.max(abs(divergence)), index)
                    gcl = J1-J0-dt*Jm*divergence
                    observe("material_GCL_relative", np.max(abs(gcl))/max(np.max(abs(J0)), 1e-30), index)
                    bulk_loss = dt*action["bulk_power"]
                    left_loss = dt*action["left_power"]
                    grad_div_work = -dt*action["grad_div_power"]
                    require(bulk_loss >= -1e-15 and left_loss >= -1e-15, "negative_physical_dissipation")
                    require(grad_div_work <= 0., "positive_independent_numerical_grad_div_work")
                    cumulative_grad_div += grad_div_work
                    observe("numerical_grad_div_dissipation_rate", action["grad_div_power"], index)
                    observe("saved_grad_div_work_discrepancy", abs(row.get("numerical_grad_div_work", 0.)-grad_div_work), index)
                    cumulative_bulk += bulk_loss; cumulative_left += left_loss
                    # Independent quadrature computation of each algebraic geometry term.
                    weights = midpoint["scalar"].W[None, :]
                    j0, j1, jm = abs(J0), abs(J1), abs(Jm)
                    ja, dj = .5*(j0+j1), j1-j0
                    mid_v = np.asarray(field)
                    change_v = np.asarray(midpoint["velocity"].interpolate(dv))
                    z0 = np.asarray(start_op["scalar"].global_coordinates())[1]
                    z1 = np.asarray(op["scalar"].global_coordinates())[1]
                    raw_wk = float(np.sum(weights*((ja-jm)*np.sum(mid_v*change_v, axis=0)
                                + .5*dj*np.sum(mid_v*mid_v, axis=0)
                                + .125*dj*np.sum(change_v*change_v, axis=0))))
                    raw_wp = float(9.81*np.sum(weights*((ja-jm)*(z1-z0)+.5*dj*(z0+z1))))
                    geometric_work = dt*action["geometric_power"]
                    hydrostatic_work = dt*action["hydrostatic_power"]
                    wk, wp = raw_wk-geometric_work, raw_wp+hydrostatic_work
                    if "raw_geometry_kinetic_work" in row:
                        for name, value in (("raw_geometry_kinetic_work", raw_wk),
                                            ("raw_geometry_potential_work", raw_wp),
                                            ("geometric_divergence_work", geometric_work),
                                            ("hydrostatic_split_work", hydrostatic_work)):
                            observe("saved_work_discrepancy", abs(row[name]-value), index)
                    cumulative_kinetic += wk; cumulative_potential += wp
                    for name, value in (("bulk_loss", bulk_loss), ("left_wall_loss", left_loss),
                                        ("geometry_kinetic_work", wk), ("geometry_potential_work", wp)):
                        observe("saved_work_discrepancy", abs(row[name]-value), index)
                    algebraic = measured["energy"]-start_energy+bulk_loss+left_loss-wk-wp-grad_div_work
                    observe("step_split_energy_relative", abs(algebraic)/E0, index)
                    # Report physical traction residuals; natural FEM traction is weak,
                    # not a pointwise Dirichlet datum to be falsely required to vanish.
                    for side in ("left", "surface"):
                        traction = _traction_norms(midpoint, vm, p, side, controls["intorder"], hydrostatic, gamma)
                        observe(side+"_traction_L2_per_density", traction["physical"], index)
                        if side == "surface":
                            observe("surface_augmented_traction_L2_per_density", traction["augmented"], index)
                            observe("surface_grad_div_traction_L2_per_density", traction["grad_div"], index)
                            observe("surface_grad_div_traction_sampled_max_per_density", traction["grad_div_sampled_max"], index)
                continuous = (measured["energy"]-E0+cumulative_bulk+cumulative_left)/E0
                split = continuous-(cumulative_kinetic+cumulative_potential+cumulative_rezone+cumulative_grad_div)/E0
                observe("continuous_energy_relative", abs(continuous), index)
                observe("split_energy_relative", abs(split), index)
                observe("cumulative_numerical_grad_div_loss_over_E0", -cumulative_grad_div/E0, index)
                observe("saved_cumulative_grad_div_work_discrepancy",
                        abs(row.get("cumulative_numerical_grad_div_work", 0.)-cumulative_grad_div), index)
                for name, value in (("cumulative_bulk_loss", cumulative_bulk),
                                    ("cumulative_left_wall_loss", cumulative_left),
                                    ("cumulative_geometry_kinetic_work", cumulative_kinetic),
                                    ("cumulative_geometry_potential_work", cumulative_potential)):
                    observe("saved_cumulative_work_discrepancy", abs(row[name]-value), index)
                observe("saved_rezone_work_discrepancy", abs(row.get("cumulative_rezone_work", 0.)-cumulative_rezone), index)
                previous_X, previous_v, previous_op = X, v, op
                previous_minimum = minimum
                previous_energy = measured["energy"]
                last_complete = {"step": index, "time_s": time_s}
                if index % 10 == 0:
                    heartbeat()
            for group in h["rejected_trials"].values():
                source = int(group.attrs["from_step"])
                require(not bool(group.attrs["accepted"]) and 0 <= source < len(keys)
                        and group.attrs["time_s"] == time_history[source], "rejected_trial_mutated_accepted_time")
            require(abs(float(h.attrs["last_accepted_time_s"])-time_history[-1]) <= 1e-12,
                    "false_final_native_time")
            require(time_history[-1] <= controls["t_end"]+1e-12, "native_exceeds_declared_horizon")
            if require_target:
                require(abs(time_history[-1]-5.) <= 1e-12, "not_real_0_to_5s_trajectory")
            from .free_boundary_geometry_audit import right_branch_gap
            resolution = []
            for index in sorted({len(keys)-1, witnesses["inverse_min_Jacobian"]}):
                geometry = h["states"][keys[index]]["geometry"][:]
                resolution.append({"step": index, "time_s": time_history[index],
                                   **right_branch_gap(geometry, triangles, boundaries, ids)})
            report.update(native_sha256=_hash_file(path), run_id=run_id,
                accepted_states=len(keys), time_start_s=time_history[0], time_end_s=time_history[-1],
                actual_5s_trajectory=abs(time_history[-1]-5.) <= 1e-12,
                full_interface_nodes=len(ids), max_curve_tessellation_bound_m=max_tessellation,
                H_final_m=H.tolist(), R_final_m=contact_history[-1], markers=previous_markers,
                split_energy_per_density={"bulk_loss": cumulative_bulk, "left_wall_loss": cumulative_left,
                    "signed_geometry_kinetic_work": cumulative_kinetic,
                    "signed_geometry_potential_work": cumulative_potential,
                    "signed_rezone_work": cumulative_rezone,
                    "signed_numerical_grad_div_work": cumulative_grad_div},
                grad_div_gamma_m2_s=gamma,
                independently_checked_rezones=rezone_count, pressure_variable=pressure_variable,
                right_branch_resolution_diagnostics=resolution,
                geometry_resolution_passed=None,
                geometry_one_contact_per_side_passed=not report["additional_solid_contact_events"],
                geometry_resolution_note="Positive exact cell Jacobians and complete curve are validity checks; full h refinement is separately required.")
            report["endpoint_divergence_note"] = "B(Xn)*vn is measured separately: the material midpoint DAE constrains B(Xmid)*vmid, not the accepted endpoint exactly. Endpoint defect requires dt refinement."
            report["grad_div_note"] = "Consistent numerical stabilization, not physical viscosity/heat. Its negative work enters only the split numerical balance; raw physical-energy acceptance remains unchanged. Physical FREE traction, augmented traction (T+gamma div(v) I)n, and the numerical contribution are separately sampled, not asserted to vanish pointwise."
            report["audit_complete"] = True
        limits = {"right_velocity_max_m_s": 1e-12, "bottom_velocity_max_m_s": 1e-12,
                  "left_normal_velocity_max_m_s": 1e-12, "P2_position_drift_m": 1e-12,
                  "solid_boundary_position_drift_m": 1e-12, "retained_running_max_error_m": 1e-12,
                  "right_H_error_m": 1e-12, "material_kinematic_defect_m": 1e-11,
                  "material_GCL_relative": 1e-9, "momentum_residual": 1e-8,
                  "weak_divergence": 1e-10, "startup_momentum": 1e-8,
                  "startup_weak_divergence": 1e-10, "startup_fixed_acceleration": 1e-12,
                  "volume_relative": 1e-6, "continuous_energy_relative": .05,
                  "split_energy_relative": 1e-6, "step_split_energy_relative": 1e-7,
                  "saved_diagnostic_discrepancy": 1e-9, "saved_work_discrepancy": 1e-10,
                  "saved_cumulative_work_discrepancy": 1e-9,
                  "saved_grad_div_work_discrepancy": 1e-10,
                  "saved_cumulative_grad_div_work_discrepancy": 1e-9,
                  "saved_rezone_work_discrepancy": 1e-9, "rezone_surface_geometry_change_m": 1e-12,
                  "rezone_constrained_velocity_m_s": 1e-12, "rezone_inverse_map_error_m": 1e-10,
                  "rezone_projection_residual": 1e-9, "rezone_momentum_error": 1e-10,
                  "rezone_weak_divergence": 1e-10, "rezone_volume_change_m2": 1e-10,
                  "rezone_saved_energy_work_discrepancy": 1e-9,
                  "rezone_saved_volume_change_discrepancy": 1e-10}
        for name, tolerance in limits.items():
            require(maxima.get(name, 0.) <= tolerance, "limit:"+name)
        report["limits"] = limits
    except Exception as error:
        # Fail closed for malformed/native-degenerate input, including third-
        # party geometry exceptions. KeyboardInterrupt/SystemExit still escape.
        require(False, "audit_exception:"+type(error).__name__+":"+str(error))
    report.update(passed=not errors and report["audit_complete"], failures=errors, maxima=maxima,
                  maximum_at_step=witnesses, last_fully_audited_state=last_complete,
                  interrupted_at_state=current_index if not report["audit_complete"] else None)
    return report
