from dataclasses import replace

import numpy as np

from sloshing.config import SimulationConfig
from sloshing.validation.qualification.comparison import compare_datasets
from sloshing.validation.qualification.dataset import run_dataset


def test_cross_mesh_comparison_and_compatible_reuse(tmp_path):
    c = SimulationConfig(nx=4, nz=6, alpha_deg=.02, integrator="sdirk2",
                         t_end=.01, dt=.00125, snapshot_dt=.005)
    m = run_dataset(c, tmp_path/"medium.h5")
    f = run_dataset(replace(c, nx=6, nz=9, mesh="fine"), tmp_path/"fine.h5")
    a = compare_datasets(m, f, tmp_path/"comparison")
    b = compare_datasets(m, f, tmp_path/"comparison")
    assert a == b
    np.testing.assert_allclose(a["integrated_domain_area"], 20., atol=1e-12)
    assert a["quadrature_is_exact_for_piecewise_fields"]
    assert a["velocity_relative_max"] > 0
    assert a["omega_full_relative_max"] > 0
    assert a["velocity_symmetry_fine_max"] < 1e-12
