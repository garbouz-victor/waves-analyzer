import numpy as np

from sloshing import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.postprocess.vorticity import VorticityProjector
from sloshing.validation.regions import VorticityRegions


def test_exact_clipped_region_integrals():
    f = FEMSystem(SimulationConfig(nx=8, nz=16))
    projector = VorticityProjector(f)
    regions = VorticityRegions(f, projector)
    x, z = projector.basis.doflocs
    constant = regions.norms(np.ones_like(x))
    np.testing.assert_allclose([constant[k]**2 for k in ("full", "wall_left", "wall_right", "surface_layer")],
                               [20, 1, 1, 2], atol=1e-12)
    linear = regions.norms(x)
    np.testing.assert_allclose(linear["wall_left"]**2, 10*(1-.9**3)/3, atol=1e-12)
    np.testing.assert_allclose(linear["wall_left"], linear["wall_right"], atol=1e-12)
