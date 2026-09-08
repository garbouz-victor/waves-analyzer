"""Exercise the actual comparison/report implementation on tiny complete histories."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from sloshing.multiphase.phase_rate_schedule import RateSchedule
from sloshing.multiphase.phase_rate_history import read_json
from sloshing.multiphase.full_rate_policy import POLICY
from test_step3a7_schedule import tiny_plan
from test_step3a8_full_phase_rate import tiny_solver
from test_step3a8_full_trajectory import make_runner


def test_full_comparison_and_figures_on_complete_tiny_data(tmp_path):
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    from sloshing.multiphase.phase_rate_trajectory import RateTrajectory
    from sloshing.multiphase.full_rate_comparison import compare
    p=read_json("validation_results/step3a4/policy/policy.json")["policy"]
    schedules=[RateSchedule(tiny_plan(),f) for f in (1,2,4)]
    for i,schedule in enumerate(schedules):
        e=IsolatedCHPhaseRate(tiny_solver())
        identity=dict(tiny_test=True,schedule_sha256=schedule.sha256)
        r=RateTrajectory(e,schedule,tmp_path/f"iso/level{i}",identity,p,float(e.metric.mass@e.phi_old.x.array),expected_crossings=0,expanded_audit=False)
        try: r.run_until(schedule.nsteps)
        finally: r.close()
    full=make_runner(tmp_path/"full")
    try:
        full.run_until(6)
        data=compare(tmp_path/"full",tmp_path/"iso",schedules,full.engine.metric.M,
            np.full(len(full.engine.phi_old.x.array),.2),-.020364675298172572,POLICY)
        assert data["no_interpolation"] and len(data["scalars"])==7
        assert len(data["fields"])==len(schedules[0].common_parent_indices)
        assert data["coupling"]["combined_D_proxy"]>=data["coupling"]["measured"]["D_coupling"]
        from sloshing.multiphase.phase_rate_history import Journal
        spec=importlib.util.spec_from_file_location("step3a8_report",Path("scripts/step3a8_report.py"))
        report=importlib.util.module_from_spec(spec); spec.loader.exec_module(report)
        report.plots(tmp_path,full.journal.rows,Journal.read(tmp_path/"iso/level0/scalar_history.jsonl"),data)
        assert len(list((tmp_path/"analysis").glob("*.pdf")))==7
    finally: full.close()
