import copy
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from test_step3a8_full_phase_rate import tiny_solver
from test_step3a7_schedule import tiny_plan
from test_step3a9_divergence import rows
from sloshing.multiphase.divergence_policy import production_policy,spatial_decision
from sloshing.multiphase.phase_rate_history import read_json,load_archive


def runner(folder,resume=False):
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    from sloshing.multiphase.divergence_trajectory import DivergenceTrajectory
    from sloshing.multiphase.phase_rate_schedule import RateSchedule
    s=tiny_solver();e=CHNSPhaseRateBE(s);schedule=RateSchedule(tiny_plan())
    p=production_policy(spatial_decision(rows()))
    return DivergenceTrajectory(e,schedule,folder,{"schedule_sha256":schedule.sha256,"test":"step3a9"},
        read_json("validation_results/step3a4/policy/policy.json")["policy"],float(e.metric.mass@s.state.x.array[e.rate_map]),
        expected_crossings=0,full_policy=p,resume=resume)


def test_new_policy_full_restart_through_dt_change(tmp_path):
    a=runner(tmp_path/"full");b=runner(tmp_path/"split")
    try:
        a.run_until(6);b.run_until(3)
        assert a.journal.rows[4]["solver"]["guess"]["sha256"]==a.journal.rows[3]["solver"]["retained_rate_sha256"]
        assert all(r["full_physical_checks"]["weak_continuity"] for r in a.journal.rows)
    finally:a.close();b.close()
    env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(Path("tests").resolve()),str(Path("src").resolve()),os.environ.get("PYTHONPATH","")]))
    code="from test_step3a9_transfer import runner; import sys; r=runner(sys.argv[1],True); r.run_until(6); r.close()"
    subprocess.run([sys.executable,"-c",code,str(tmp_path/"split")],env=env,check=True,timeout=180)
    x,mx=load_archive(tmp_path/"full/checkpoints/000006");y,my=load_archive(tmp_path/"split/checkpoints/000006")
    for key in ("state","old","older","phi","mu","u","pi","phase_rate"):
        assert np.allclose(x[key],y[key],rtol=1e-10,atol=1e-13)
    for key in ("CH","viscous","slip"):
        assert mx["diagnostics"]["cumulative"][key]==pytest.approx(my["diagnostics"]["cumulative"][key],rel=1e-10,abs=1e-18)


@pytest.mark.parametrize("failure",["weak","nonfinite","SNES","material"])
def test_new_trajectory_constraint_failure_transaction(tmp_path,monkeypatch,failure):
    from sloshing.multiphase.phase_rate_trajectory import ScientificFailure
    r=runner(tmp_path)
    try:
        r.run_until(2);before=r.physical_arrays();cum=copy.deepcopy(r.diag.cumulative)
        checkpoint=(tmp_path/"LATEST.json").read_bytes();history=(tmp_path/"scalar_history.jsonl").read_bytes()
        if failure in ("weak","nonfinite"):
            measure=r.full_observer.projection.measure
            def bad(weak):
                row=measure(weak)
                row["weak_continuity_Riesz" if failure=="weak" else "strong_divergence_L2"]=2e-12 if failure=="weak" else np.nan
                return row
            monkeypatch.setattr(r.full_observer.projection,"measure",bad)
        else:
            def fail(*a,**kw):raise RuntimeError("injected "+failure)
            if failure=="SNES":monkeypatch.setattr(r.engine.problem,"solve",fail)
            else:
                import sloshing.multiphase.nonlinear_accuracy as n
                monkeypatch.setattr(n,"check_candidate_material",fail)
        with pytest.raises(ScientificFailure):r.run_until(6)
        assert r.s.step_number==2 and cum==r.diag.cumulative
        assert all(np.array_equal(v,r.physical_arrays()[k]) for k,v in before.items())
        assert checkpoint==(tmp_path/"LATEST.json").read_bytes() and history==(tmp_path/"scalar_history.jsonl").read_bytes()
        assert not (tmp_path/"COMPLETE.json").exists()
    finally:r.close()


def test_production_audit_dependency_and_schedule_identity():
    from sloshing.multiphase.benchmarks import divergence_transfer as b
    assert b.audit_proof()["decision"]["passed"]
    assert b.schedules()[0].sha256==b.L0_SHA
    assert b.POLICY["version"]=="step3a9-v1" and "strong_divergence_max" not in b.POLICY
