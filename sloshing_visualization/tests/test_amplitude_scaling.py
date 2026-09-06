from dataclasses import replace

import numpy as np

from sloshing.config import SimulationConfig
from sloshing.solver import SloshingSolver


def test_real_production_path_scales_with_tan_alpha():
    base = SimulationConfig(nx=8, nz=16, alpha_deg=.2, integrator="sdirk2",
                            t_end=.025, dt=.00125, snapshot_dt=.005)
    small = replace(base, alpha_deg=.02)
    a, b = SloshingSolver(base), SloshingSolver(small)
    sa, sb = a.integrator.initial_state(), b.integrator.initial_state()
    scale = small.slope/base.slope
    for step in range(base.nsteps+1):
        if step:
            sa, sb = a.integrator.advance(sa), b.integrator.advance(sb)
        # Compare fields, not elementwise ratios which divide by physical zeros.
        np.testing.assert_allclose(sb.velocity, scale*sa.velocity, atol=2e-14, rtol=2e-10)
        np.testing.assert_allclose(sb.eta, scale*sa.eta, atol=2e-15, rtol=2e-10)
        np.testing.assert_allclose(b.fem.energy(sb.velocity, sb.eta),
                                   scale**2*np.array(a.fem.energy(sa.velocity, sa.eta)),
                                   atol=1e-22, rtol=2e-10)
