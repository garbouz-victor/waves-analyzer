import copy
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from sloshing.multiphase.phase_rate_history import load_archive,read_json
from sloshing.multiphase.phase_rate_schedule import RateSchedule
from sloshing.multiphase.full_rate_policy import POLICY
from test_step3a7_schedule import tiny_plan
from test_step3a8_full_phase_rate import tiny_solver


def make_runner(folder,resume=False):
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    from sloshing.multiphase.full_rate_trajectory import FullRateTrajectory
    s=tiny_solver(); e=CHNSPhaseRateBE(s)
    schedule=RateSchedule(tiny_plan()); identity=dict(tiny_test=True,schedule_sha256=schedule.sha256)
    p=read_json("validation_results/step3a4/policy/policy.json")["policy"]
    mass=float(e.metric.mass@s.state.x.array[e.rate_map])
    return FullRateTrajectory(e,schedule,folder,identity,p,mass,expected_crossings=0,
        resume=resume,full_policy=POLICY)


def test_coupled_multistep_changing_dt_and_restart(tmp_path):
    a=make_runner(tmp_path/"full"); b=make_runner(tmp_path/"split")
    try:
        result=a.run_until(6); b.run_until(3)
        assert result["execution_status"]=="complete"
        rows=a.journal.rows
        assert "FULL" in rows[0]["benchmark"]
        assert rows[4]["solver"]["guess"]["sha256"]==rows[3]["solver"]["retained_rate_sha256"]
        assert rows[4]["dt"]==2*rows[3]["dt"]
    finally: a.close(); b.close()
    code="from test_step3a8_full_trajectory import make_runner; import sys; r=make_runner(sys.argv[1],True); r.run_until(6); r.close()"
    env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(Path("tests").resolve()),str(Path("src").resolve()),os.environ.get("PYTHONPATH","")]))
    subprocess.run([sys.executable,"-c",code,str(tmp_path/"split")],check=True,env=env,timeout=180)
    x,mx=load_archive(tmp_path/"full/checkpoints/000006"); y,my=load_archive(tmp_path/"split/checkpoints/000006")
    for k in ("u","pi","phi","mu","phase_rate","state","old","older"):
        assert np.allclose(x[k],y[k],atol=1e-13,rtol=1e-10),k
    for k in ("CH","viscous","slip"):
        assert mx["diagnostics"]["cumulative"][k]==pytest.approx(my["diagnostics"]["cumulative"][k],abs=1e-18,rel=1e-10)


@pytest.mark.parametrize("kind",["SNES","material","physical"])
def test_coupled_failed_candidate_keeps_state_history_and_checkpoint(tmp_path,monkeypatch,kind):
    from sloshing.multiphase.phase_rate_trajectory import ScientificFailure
    r=make_runner(tmp_path)
    try:
        r.run_until(2); arrays=r.physical_arrays(); diag=copy.deepcopy(r.diag.cumulative)
        attrs={k:getattr(r.s,k) for k in ("time","step_number","current_dt","current_phase")}
        journal=(tmp_path/"scalar_history.jsonl").read_bytes(); checkpoint=(tmp_path/"LATEST.json").read_bytes()
        def fail(*args,**kwargs): raise RuntimeError("Injected "+kind)
        if kind=="SNES": monkeypatch.setattr(r.engine.problem,"solve",fail)
        elif kind=="material":
            import sloshing.multiphase.nonlinear_accuracy as h
            monkeypatch.setattr(h,"check_candidate_material",fail)
        else: monkeypatch.setattr(r.full_observer,"measure",fail)
        with pytest.raises(ScientificFailure): r.run_until(6)
        assert all(np.array_equal(v,r.physical_arrays()[k]) for k,v in arrays.items())
        assert all(getattr(r.s,k)==v for k,v in attrs.items()) and r.diag.cumulative==diag
        assert journal==(tmp_path/"scalar_history.jsonl").read_bytes()
        assert checkpoint==(tmp_path/"LATEST.json").read_bytes()
        assert not (tmp_path/"COMPLETE.json").exists()
    finally: r.close()
