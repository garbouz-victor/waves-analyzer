"""Cheap container-only coupled tests; no production tank in ordinary CI."""
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip("dolfinx",reason="STEP 3 coupled tests run in pinned FEniCSx container")
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.initialization import flat_interface, sessile_drop
from sloshing.multiphase.diagnostics import Diagnostics
from sloshing.multiphase.storage import write_checkpoint, load_checkpoint
from sloshing.multiphase.interface import contour_segments, connected_polylines, fit_circle


def tiny_config():
    path=Path(__file__).resolve().parents[1]/"configs/step3/benchmark_flat.json"
    return ModelConfig.from_json(path).changed(nx=8,nz=8,t_end=.004)


def test_step3_flat_coupled_mass_and_energy():
    c=tiny_config()
    solver=CHNSSolver(c)
    solver.initialize(flat_interface(c.epsilon))
    d=Diagnostics(solver)
    first=d.measure()
    for _ in range(2):
        solver.advance();last=d.measure()
    assert abs(last["mass_error_relative_to_domain"])<1e-9
    assert last["E_total"] <= first["E_total"]+1e-10
    assert last["wall_normal_velocity_L2"]<1e-14
    assert solver.logs[-1]["residual"]<1e-8


def test_step3_checkpoint_preserves_bdf_history(tmp_path):
    c=tiny_config()
    solver=CHNSSolver(c);solver.initialize(flat_interface(c.epsilon))
    diag=Diagnostics(solver);diag.measure()
    solver.advance();diag.measure()
    write_checkpoint(tmp_path,solver,diag,status="running")
    solver.advance();uninterrupted=diag.measure()
    expected=solver.state.x.array.copy()
    restored=CHNSSolver(c);drestored=Diagnostics(restored)
    with pytest.raises(ValueError,match="not complete"):
        load_checkpoint(tmp_path,restored,drestored)
    load_checkpoint(tmp_path,restored,drestored,allow_running=True)
    restored.advance();result=drestored.measure()
    assert np.max(np.abs(expected-restored.state.x.array))<1e-10
    assert result["energy_budget_defect"]==pytest.approx(uninterrupted["energy_budget_defect"],abs=1e-12)


@pytest.mark.parametrize("angle",[60.,120.])
def test_step3_wetting_coupled_smoke(angle):
    c=tiny_config().changed(z_min=0,z_max=1,theta_equilibrium_deg=angle,
        no_slip_walls=(),free_slip_walls=("top",),wetting_walls=("bottom",))
    solver=CHNSSolver(c)
    solver.initialize(sessile_drop(c.epsilon,.3,90))
    d=Diagnostics(solver);first=d.measure()
    solver.advance();last=d.measure()
    assert abs(last["mass_error_relative_to_domain"])<1e-9
    assert last["E_total"] < first["E_total"]
    # This is a nonlinear/wetting SMOKE test, not an angle-convergence claim.


def test_step3_fem_contour_exact_quadratic():
    c=tiny_config()
    solver=CHNSSolver(c)
    # Not a PDE benchmark: a P2 polynomial tests actual FE edge evaluation.
    solver.state.sub(2).interpolate(lambda x: .3**2-x[0]**2-x[1]**2)
    solver.state.x.scatter_forward()
    segments=contour_segments(solver)
    assert len(connected_polylines(segments))==1
    fitted=fit_circle(segments)
    assert fitted["radius"]==pytest.approx(.3,abs=1e-10)
    assert fitted["radial_rms"]<1e-10


def test_step3_sessile_contour_at_mesh_vertex():
    path=Path(__file__).resolve().parents[1]/"configs/step3/benchmark_contact.json"
    c=ModelConfig.from_json(path).changed(nx=16,nz=8,refinement_levels=0,epsilon=.08)
    solver=CHNSSolver(c)
    solver.state.sub(2).interpolate(sessile_drop(c.epsilon,.3,90.,0.))
    solver.state.x.scatter_forward()
    assert len(connected_polylines(contour_segments(solver)))==1
