"""Exact P2 slopes and regional FEM derivative diagnostics."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from ...diagnostics import surface_slopes
from ..contact import surface_value


def slope_regions(x, eta):
    left, right = surface_slopes(x, eta)
    edges = x[::2]
    def interval(lo, hi):
        lo = np.maximum(edges[:-1], lo)
        hi = np.minimum(edges[1:], hi)
        hit = hi >= lo
        if not np.any(hit):
            return 0.
        derivative = (right-left)/np.diff(edges)
        dl = left+derivative*(lo-edges[:-1])
        dr = left+derivative*(hi-edges[:-1])
        return float(max(np.max(abs(dl[hit])), np.max(abs(dr[hit]))))
    a = max(abs(x[0]), abs(x[-1]))
    return {"slope_global": interval(-a, a), "slope_bulk": interval(-.9*a, .9*a),
            "slope_intermediate": interval(-.98*a, .98*a),
            "slope_contact": max(interval(-a, -.98*a), interval(.98*a, a))}


def modal_projection(x, eta, a):
    nodes, weights = np.polynomial.legendre.leggauss(8)
    edges = x[::2]
    half = np.diff(edges)/2
    xx = (edges[:-1]+edges[1:])[:, None]/2+half[:, None]*nodes
    phi = np.sin(np.pi*xx/(2*a))
    # ||sin(k1 x)||^2 on [-a,a] equals a; no inviscid dynamics imposed.
    return float(np.sum(half[:, None]*weights*surface_value(x, eta, xx)*phi)/a)


def pinned_discrete_equilibrium(x, eta_endpoints, g):
    """Algebraic interpretation ONLY, not an independent physical benchmark.

    For antisymmetric pinned endpoints, q=0 and v=0 admit S_interior*eta=0.
    The resulting localized P2 trace is in ker(C); its gradient is mesh-based.
    It is not a continuum wall film and need not be reached by t=5.
    """
    if not np.isclose(sum(eta_endpoints), 0., atol=1e-15):
        raise ValueError("This diagnostic assumes antisymmetric zero-mean release")
    h = np.diff(x[::2])
    ids = 2*np.arange(len(h))[:, None]+np.arange(3)
    local = h[:, None, None]/30*np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]])
    S = coo_matrix((local.ravel(), (np.broadcast_to(ids[:, :, None], local.shape).ravel(),
                                   np.broadcast_to(ids[:, None, :], local.shape).ravel())),
                   shape=(len(x), len(x))).tocsr()
    eta = np.zeros_like(x)
    eta[[0, -1]] = eta_endpoints
    eta[1:-1] = spsolve(S[1:-1, 1:-1], -(S@eta)[1:-1])
    return eta, {"purpose": "discrete algebraic interpretation; not an independent benchmark or a wall-film model",
                 "max_slope": slope_regions(x, eta)["slope_global"],
                 "potential_energy": float(g/2*eta@(S@eta)),
                 "interior_surface_load_residual": float(np.max(abs((S@eta)[1:-1]))),
                 "surface_integral": float(np.ones(len(x))@(S@eta))}


class DerivativeMetrics:
    def __init__(self, fem, regions):
        self.fem = fem
        self.bulk_mass = regions.matrices["full"]-regions.matrices["wall_left"]-regions.matrices["wall_right"]
        coords = fem.velocity.global_coordinates()
        vertices = fem.mesh.p[:, fem.mesh.t].transpose(2, 1, 0)
        # Each P2 derivative is affine. Recover its vertex values from exact
        # FEM quadrature samples (no finite differences or cross-element average).
        design = np.stack((np.ones_like(coords[0]), coords[0], coords[1]), axis=-1)
        self.vertex_map = np.einsum("evi,eiq->evq", np.stack((np.ones_like(vertices[:, :, 0]), vertices[:, :, 0], vertices[:, :, 1]), axis=-1),
                                    np.linalg.pinv(design))

    def evaluate(self, full):
        f = self.fem
        grad = f.velocity.interpolate(full).grad
        div = grad[0, 0]+grad[1, 1]
        values = np.einsum("evq,eq->ev", self.vertex_map, div)
        div_bulk = np.sqrt(max(0., np.einsum("ei,eij,ej->", values, self.bulk_mass, values)))
        gradient = np.einsum("evq,abeq->abev", self.vertex_map, grad)
        grad_bulk = np.sqrt(max(0., np.einsum("abei,eij,abej->", gradient, self.bulk_mass, gradient)))
        return {"divergence_bulk_l2": float(div_bulk), "gradient_bulk_l2": float(grad_bulk),
                "divergence_bulk_relative": float(div_bulk/max(grad_bulk, 1e-30))}
