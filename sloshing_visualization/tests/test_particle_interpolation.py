import numpy as np
import pytest

from sloshing.config import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.postprocess.vorticity import VorticityProjector
from sloshing.validation.qualification.fem_evaluation import FEMGeometry
from sloshing.validation.qualification.overlay import common_quadrature


@pytest.mark.parametrize("nx", [4, 5])
def test_actual_p2_interpolation_matches_fem_probes(nx):
    f = FEMSystem(SimulationConfig(nx=nx, nz=8))
    geom = FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs)
    rng = np.random.default_rng(431)
    points = np.vstack((rng.uniform(-1, 1, 100), rng.uniform(-10, 0, 100)))
    coeff = rng.normal(size=f.scalar.N)
    np.testing.assert_allclose(geom.evaluate(coeff, points), f.scalar.probes(points)@coeff, atol=5e-12)


def test_overlay_integrates_both_p2_fields_without_visualization_sampling():
    aa = FEMSystem(SimulationConfig(a=1, d=1, nx=4, nz=4, x_grading=1.1, z_grading=.7))
    bb = FEMSystem(SimulationConfig(a=1, d=1, nx=6, nz=6, x_grading=1.1, z_grading=.7))
    a, b = [FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs) for f in (aa, bb)]
    points, weights, _ = common_quadrature(a, b)
    np.testing.assert_allclose(weights.sum(), 2., atol=1e-13)
    exact = lambda x, z: x*x+2*z*z+x*z
    av = a.evaluate(exact(*aa.scalar.doflocs), points)
    bv = b.evaluate(exact(*bb.scalar.doflocs), points)
    np.testing.assert_allclose(av, bv, atol=5e-13)
    # Independently integrate polynomial squared by tensor Gauss quadrature.
    q, w = np.polynomial.legendre.leggauss(5)
    reference = np.sum(w[:, None]*w[None, :]/2*exact(q[:, None], (q[None, :]-1)/2)**2)
    np.testing.assert_allclose(np.sum(weights*av**2), reference, atol=1e-12)


def test_overlay_resolves_discontinuous_vorticity_across_both_meshes():
    aa = FEMSystem(SimulationConfig(a=1, d=1, nx=4, nz=4))
    bb = FEMSystem(SimulationConfig(a=1, d=1, nx=6, nz=6))
    projectors = [VorticityProjector(f) for f in (aa, bb)]
    geoms = [FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs, p.basis.element_dofs) for f, p in zip((aa, bb), projectors)]
    rng = np.random.default_rng(54)
    fields = [rng.normal(size=p.basis.N) for p in projectors]
    norms = []
    for order in (4, 6):
        points, weights, _ = common_quadrature(*geoms, order=order)
        values = [g.evaluate(v, points, 1) for g, v in zip(geoms, fields)]
        norms.append(np.sum(weights*(values[0]-values[1])**2))
        for value, coeff, p in zip(values, fields, projectors):
            reference = np.sum(p.basis.interpolate(coeff)**2*p.basis.dx)
            np.testing.assert_allclose(np.sum(weights*value**2), reference, atol=2e-12)
    np.testing.assert_allclose(*norms, atol=3e-12)
