"""Tiny meshes only. Production trajectories are never part of CI."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from sloshing.multiphase.phase_rate_history import read_json,load_archive,Journal
from sloshing.multiphase.phase_rate_schedule import RateSchedule
from test_step3a7_schedule import tiny_plan


def make_runner(folder,factor=1,resume=False,dt=1e-6):
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    from sloshing.multiphase.phase_rate_trajectory import RateTrajectory
    c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4,nz=4,theta_equilibrium_deg=90.,g=0.,a_x=0.,epsilon=1.,quadrature_degree=12)
    s=CHNSSolver(c)
    engine=IsolatedCHPhaseRate(s,lambda x:.2+.001*np.cos(2*x[0])*np.cos(3*x[1]))
    mass=float(engine.metric.mass @ s.state.x.array[engine.physical_phi_map])
    policy=read_json("validation_results/step3a4/policy/policy.json")["policy"]
    schedule=RateSchedule(tiny_plan(dt),factor)
    identity={"tiny_test_only":True,"schedule_sha256":schedule.sha256,"source":"test-fixed"}
    return RateTrajectory(engine,schedule,folder,identity,policy,mass,
        expected_crossings=0,resume=resume,expanded_audit=True)


def test_multistep_block_transition_and_retained_unscaled_guess(tmp_path):
    runner=make_runner(tmp_path)
    try:
        result=runner.run_until(6)
        assert result["execution_status"]=="complete"
        rows=runner.journal.rows
        for before,after in zip(rows[1:],rows[2:]):
            assert after["solver"]["guess"]["sha256"]==before["solver"]["retained_rate_sha256"]
        assert rows[4]["dt"]==2*rows[3]["dt"]
        assert all(all(r["physical_gates"].values()) for r in rows)
        assert result["min_certified_cells"]>=8
    finally: runner.close()


def test_restart_in_new_process_across_dt_boundary(tmp_path):
    full=make_runner(tmp_path/"full"); split=make_runner(tmp_path/"split")
    try:
        full.run_until(6); split.run_until(3)
    finally: full.close(); split.close()
    # Fresh Python/FEniCSx runner; the resumed first interval changes dt.
    code="from test_step3a7_trajectory import make_runner; import sys; r=make_runner(sys.argv[1],resume=True); r.run_until(6); r.close()"
    env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(Path("tests").resolve()),str(Path("src").resolve()),os.environ.get("PYTHONPATH","")]))
    subprocess.run([sys.executable,"-c",code,str(tmp_path/"split")],check=True,env=env,timeout=180)
    a,ma=load_archive(tmp_path/"full/checkpoints/000006")
    b,mb=load_archive(tmp_path/"split/checkpoints/000006")
    for key in ("phi","mu","phase_rate","state","old","older"):
        assert np.allclose(a[key],b[key],atol=1e-13,rtol=1e-12)
    assert ma["time"]==mb["time"] and ma["accepted_steps"]==mb["accepted_steps"]==6
    assert ma["diagnostics"]["cumulative"]["CH"]==pytest.approx(mb["diagnostics"]["cumulative"]["CH"],rel=1e-12)
    for left,right in zip(Journal.read(tmp_path/"full/scalar_history.jsonl"),Journal.read(tmp_path/"split/scalar_history.jsonl")):
        for key in ("time","E_total","CH_dissipation","cumulative_CH_dissipation","energy_budget_defect"):
            assert left[key]==pytest.approx(right[key],abs=1e-14,rel=1e-12)


@pytest.mark.parametrize("kind",["SNES","material","work"])
def test_failed_step_preserves_accepted_history_rate_and_checkpoint(tmp_path,monkeypatch,kind):
    from sloshing.multiphase.phase_rate_trajectory import ScientificFailure
    runner=make_runner(tmp_path)
    try:
        runner.run_until(2)
        before=runner.physical_arrays(); attrs={k:getattr(runner.s,k) for k in ("time","step_number","current_dt","current_phase")}
        diag=copy.deepcopy(runner.diag.cumulative); previous=copy.deepcopy(runner.diag.previous)
        history=(tmp_path/"scalar_history.jsonl").read_bytes(); latest=(tmp_path/"LATEST.json").read_bytes()
        def fail(*args,**kwargs): raise RuntimeError("injected "+kind)
        if kind=="SNES": monkeypatch.setattr(runner.engine.problem,"solve",fail)
        elif kind=="material":
            import sloshing.multiphase.nonlinear_accuracy as h
            monkeypatch.setattr(h,"check_candidate_material",fail)
        else: monkeypatch.setattr(runner.work,"measure",fail)
        with pytest.raises(ScientificFailure): runner.run_until(6)
        assert all(np.array_equal(a,runner.physical_arrays()[k]) for k,a in before.items())
        assert all(getattr(runner.s,k)==v for k,v in attrs.items())
        assert runner.diag.cumulative==diag and runner.diag.previous==previous
        assert np.array_equal(runner.observer.state.x.array,runner.s.state.x.array)
        assert (tmp_path/"scalar_history.jsonl").read_bytes()==history
        assert (tmp_path/"LATEST.json").read_bytes()==latest
        assert not (tmp_path/"COMPLETE.json").exists()
    finally: runner.close()
    with pytest.raises(ValueError,match="failure"): make_runner(tmp_path,resume=True)


def test_tiny_three_level_nonlinear_convergence(tmp_path):
    from sloshing.multiphase.phase_rate_schedule import scalar_convergence
    I=[]; F=[]; B=[]; energy=[]; fields=[]
    schedules=[]
    for i,factor in enumerate((1,2,4)):
        runner=make_runner(tmp_path/f"level{i}",factor,dt=1e-4)
        try:
            status=runner.run_until(runner.schedule.nsteps)
            I.append(status["integrated_D_CH"]); F.append(status["DeltaF"])
            B.append(status["final_budget_defect"]); energy.append(status["final_F"])
            fields.append(runner.physical_arrays()["phi"])
            metric=runner.engine.metric; policy=runner.policy; schedules.append(runner.schedule)
        finally: runner.close()
    errors=[metric.norm(fields[k+1]-fields[k]) for k in (0,1)]
    result=scalar_convergence(I,F,B,errors,energy,policy)
    assert result["gates"]["budget_decreases"]
    assert result["gates"]["endpoint_phi_converges"]
    assert result["gates"]["dissipation_gaps_decrease"]
    from sloshing.multiphase.phase_rate_analysis import compare_complete_levels
    analysis=compare_complete_levels(tmp_path,schedules,metric.M,policy)
    assert np.allclose(analysis["endpoint_phi_gaps"],errors)
    assert len(analysis["common_times"])==len(schedules[0].common_parent_indices)
    assert analysis["levels"][2]["accepted_steps"]==24
    # Exercise the complete analysis/figure implementation before source freeze.
    import importlib.util
    spec=importlib.util.spec_from_file_location("step3a7_report",Path("scripts/step3a7_report.py"))
    report=importlib.util.module_from_spec(spec); spec.loader.exec_module(report)
    report.plots(tmp_path,schedules,analysis,{})
    assert (tmp_path/"analysis/energy.pdf").is_file()
