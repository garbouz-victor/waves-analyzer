import h5py
import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver
from sloshing.postprocess.vorticity import VorticityProjector
from sloshing.storage import SnapshotWriter, save_run


def test_vorticity_sign_and_exact_fem_derivative(short_run):
    solver, _ = short_run
    f = solver.fem
    full = np.zeros(f.velocity.N)
    x, z = f.scalar.doflocs
    full[f.component_dofs[0]] = x*z
    full[f.component_dofs[1]] = x*x-z*z
    projector = VorticityProjector(f)
    omega = projector.evaluate(full)
    np.testing.assert_allclose(omega, projector.basis.doflocs[0], atol=5e-11)


def test_hdf5_roundtrip_and_csv(tmp_path):
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, t_end=.025,
                                           visualization_nx=9, visualization_nz=13))
    path = save_run(solver, tmp_path / "sample.h5")
    with h5py.File(path) as h:
        assert h.attrs["status"] == "complete"
        np.testing.assert_array_equal(h["times"][:], [0, .025])
        assert h["fields/u"].shape == (2, solver.fem.scalar.N)
        assert h["fields/omega"].shape == (2, 3*solver.fem.mesh.nelements)
        assert h["visualization/u"].shape == (2, 13, 9)
        assert not np.any(h["fields/u"][0])
        assert np.max(abs(h["fields/q"][0])) > .01  # Initial non-equilibrium pressure is solved.
        for name in ("u", "w"):
            np.testing.assert_allclose(h[f"visualization/{name}"][:, :, [0, -1]], 0, atol=1e-12)
    lines = (tmp_path / "sample" / "diagnostics.csv").read_text().splitlines()
    assert len(lines) == 3
    assert "left_w_max" in lines[0]
    with pytest.raises(FileExistsError):
        save_run(solver, path)


def test_failed_write_is_marked_and_closed(tmp_path, short_run):
    solver, snapshots = short_run
    path = tmp_path / "failed.h5"
    with pytest.raises(RuntimeError, match="intentional failure"):
        with SnapshotWriter(path, solver, visualization_grid=False) as writer:
            writer.append(*snapshots[0])
            raise RuntimeError("intentional failure")
    with h5py.File(path) as h:
        assert h.attrs["status"] == "failed"
        assert "intentional failure" in h.attrs["failure"]
