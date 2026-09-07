"""Tiny repairs smoke/restart/actual-FEM resolution tests; not full benchmarks."""
from pathlib import Path
import numpy as np
import pytest
pytest.importorskip("dolfinx",reason="STEP 3A.1 FEM tests require the pinned container")
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.diagnostics import Diagnostics,interface_resolution
from sloshing.multiphase.initialization import flat_interface
from sloshing.multiphase.storage import write_checkpoint,load_checkpoint


def config(**changes):
    path=Path(__file__).resolve().parents[1]/"configs/step3/benchmark_flat.json"
    return ModelConfig.from_json(path).changed(nx=6,nz=6,**changes)


@pytest.mark.parametrize("checkpoint_step",[1,3,4,5])
def test_restart_during_startup_switch_and_bdf2(tmp_path,checkpoint_step):
    c=config(startup_dt=.0005,startup_t_end=.0015,dt=.002,t_end=.0075)
    solver=CHNSSolver(c);solver.initialize(flat_interface(c.epsilon))
    diag=Diagnostics(solver);diag.measure()
    for _ in range(checkpoint_step):
        solver.advance();diag.measure()
    write_checkpoint(tmp_path,solver,diag,"running")
    clone=CHNSSolver(c);clone_diag=Diagnostics(clone)
    load_checkpoint(tmp_path,clone,clone_diag,allow_running=True)
    assert clone.current_phase==solver.current_phase
    assert clone.history_spacing==solver.history_spacing
    for _ in range(checkpoint_step,solver.schedule.nsteps):
        solver.advance();expected=diag.measure()
        clone.advance();actual=clone_diag.measure()
    assert clone.state.x.array==pytest.approx(solver.state.x.array,abs=1e-10)
    assert actual["energy_budget_defect"]==pytest.approx(expected["energy_budget_defect"],abs=1e-12)
    assert abs(actual["mass_error_relative_to_domain"])<1e-9


@pytest.mark.parametrize("degree",[1,2])
def test_actual_fem_band_crossing_near_corner_not_lost(degree):
    solver=CHNSSolver(config(phase_degree=degree))
    solver.state.sub(2).interpolate(lambda x:100*(x[0]+.5)+100*(x[1]+.5)-1.)
    solver.state.x.scatter_forward()
    measured=interface_resolution(solver)
    assert measured["interface_present"] and measured["active_cell_count"]>0
    assert measured["worst_cell"]["bounds"][0]==pytest.approx([-.5,-.5])


def test_box_refinement_improves_actual_fem_resolution():
    c=config(refinement_boxes=({"x_min":-.5,"x_max":.5,"z_min":-.2,"z_max":.2,"levels":1},))
    coarse=CHNSSolver(config());fine=CHNSSolver(c)
    counts=[]
    for solver in (coarse,fine):
        solver.state.sub(2).interpolate(flat_interface(solver.config.epsilon))
        solver.state.x.scatter_forward()
        counts.append(interface_resolution(solver)["cells_across_transition_min"])
    assert counts[1]>1.5*counts[0]


def test_failed_run_retains_every_accepted_power_sample(tmp_path,monkeypatch):
    import json
    from sloshing.multiphase.benchmarks.common import run_case
    from sloshing.multiphase.energy_validation import read_history
    original=CHNSSolver.advance
    def injected_failure(solver):
        if solver.step_number==1:
            raise RuntimeError("injected failure after one accepted interval")
        return original(solver)
    monkeypatch.setattr(CHNSSolver,"advance",injected_failure)
    with pytest.raises(RuntimeError,match="injected"):
        run_case(config(),flat_interface(.05),tmp_path)
    rows=read_history(tmp_path/"history.csv")
    assert len(rows)==2 and rows[0]["time"]==0.
    assert "CH_dissipation" in rows[0] and "local_budget_defect" in rows[1]
    assert json.loads((tmp_path/"summary.json").read_text())["status"]=="failed"


def test_historical_outputs_cannot_be_overwritten():
    from sloshing.multiphase.benchmarks.common import run_case
    root=Path(__file__).resolve().parents[1]
    with pytest.raises(ValueError,match="read-only"):
        run_case(config(),flat_interface(.05),root/"validation_results/step3/forbidden_new_run")
