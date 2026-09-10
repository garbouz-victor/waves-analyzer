"""Bounded tests of the declared law and exact checkpoint continuation."""
import json

import h5py
import numpy as np
import pytest

from sloshing.pinned_wetting.navier import NavierConfig, NavierFEM, NavierIntegrator, channel_benchmark
from sloshing.pinned_wetting.extension_run import execute_case


def test_navier_zero_disturbance_is_exact_equilibrium():
    c = NavierConfig(nx=8,nz=12,alpha_deg=0.,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2")
    integrator = NavierIntegrator(NavierFEM(c))
    s = integrator.initial_state()
    for _ in range(c.nsteps):
        s = integrator.advance(s)
    assert np.max(abs(s.velocity)) == np.max(abs(s.eta)) == 0.
    assert s.wall_dissipated_energy == s.bulk_dissipated_energy == s.rk_energy_correction == 0.


@pytest.mark.parametrize("b", [.25,.5,1.])
def test_dimensionful_slip_channel_benchmark(b):
    c = NavierConfig(nx=8,nz=12,slip_length_m=b,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2")
    assert channel_benchmark(c)["relative_error"] < 1e-10


def test_free_endpoint_and_impermeable_sides_no_slip_bottom():
    c = NavierConfig(nx=8,nz=12,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2")
    f = NavierFEM(c)
    assert [f.R.getrow(i).nnz for i in (0,f.R.shape[0]-1)] == [1,1]
    integrator = NavierIntegrator(f)
    initial = integrator.initial_state()
    advanced = integrator.advance(initial)
    full = f.expand_velocity(advanced.velocity)
    assert np.array_equal(full[f.wall_dofs["bottom"]], np.zeros(len(f.wall_dofs["bottom"])))
    for side in ("left","right"):
        normals = np.intersect1d(f.wall_dofs[side],f.component_dofs[0])
        assert np.array_equal(full[normals],np.zeros(len(normals)))
    assert np.min(abs(advanced.eta[[0,-1]]-initial.eta[[0,-1]])) > 1e-9
    assert advanced.wall_dissipated_energy > 0 and advanced.bulk_dissipated_energy > 0


def test_real_navier_resume_preserves_created_record_and_all_stage_fields(tmp_path):
    c = NavierConfig(nx=8,nz=12,dt=.01,t_end=1.,snapshot_dt=.01,integrator="sdirk2")
    identity = {"run_id":"TEST_REAL_NAVIER_CONTINUATION","source_hash":"same_test_source"}
    first, split = tmp_path/"uninterrupted", tmp_path/"resumed"
    first.mkdir(); split.mkdir()
    noop = lambda *args: None
    execute_case(tmp_path, first, c, identity, {}, noop)
    partial = execute_case(tmp_path, split, c, identity, {}, noop, stop_after=85)
    assert not partial["completed"] and partial["accepted_step"] == 85
    with h5py.File(split/"native.h5") as h:
        assert len(json.loads(h["accepted/85"].attrs["coating"])["markers"]) >= 3
    execute_case(tmp_path, split, c, identity, {}, noop)
    with h5py.File(first/"native.h5") as a, h5py.File(split/"native.h5") as b:
        for key in a["accepted"]:
            aa, bb = a["accepted"][key], b["accepted"][key]
            assert set(aa) == set(bb)
            for field in aa:
                np.testing.assert_array_equal(aa[field][:],bb[field][:])
            for field in ("coating","integrator_scalars","state_id"):
                assert aa.attrs[field] == bb.attrs[field]
    with pytest.raises(ValueError,match="identity"):
        execute_case(tmp_path,split,c,{**identity,"source_hash":"changed"},{},noop)
