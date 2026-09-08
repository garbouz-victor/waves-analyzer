import numpy as np
import pytest
from sloshing.multiphase.coupling_transfer import scalar_coupling,field_coupling,transfer_gates
from sloshing.multiphase.full_rate_policy import POLICY,L0_SHA,cost_forecast,session_result


def synthetic(cpl=1e-6):
    scalar=scalar_coupling(dict(D=1.+cpl,F=1.+cpl,kinetic=1e-8,visc=1e-8,slip=1e-8,mass=0.),
        dict(D=1.,F=1.,mass=0.),1.,1.,POLICY)
    field=field_coupling(cpl,1e-5,1e-5+cpl,.1,1.,cpl,.1,1.,POLICY)
    final=dict(I_hydro=1e-8,excess=1.,I_CH_full=1.+cpl,I_iso0=1.,I_iso2=1.00001,
        F_CH_full=1.+cpl,F_iso0=1.,F_iso2=1.00001,budget=1e-3,DeltaE=-1.)
    return scalar,field,final


def test_coupling_less_than_uncertainty_and_triangle_proxy():
    s,f,last=synthetic(); r=transfer_gates([s],[f],last,POLICY)
    assert r["passed"]
    assert r["combined_D_proxy"]==pytest.approx(1.1e-5)
    assert not r["rigorous_full_continuum_bound"]
    s,f,last=synthetic(2e-5); r=transfer_gates([s],[f],last,POLICY)
    assert not r["passed"] and not r["gates"]["phi_coupling"]
    assert not r["gates"]["D_coupling"] and not r["gates"]["F_coupling"]


def test_significance_floors_and_nonfinite():
    s=scalar_coupling(dict(D=1e-20,F=1.,kinetic=0.,visc=1e-25,slip=0.,mass=0.),
        dict(D=1e-21,F=1.,mass=0.),1.,1.,POLICY)
    assert s["CH_power_relative_significant"] is None and s["hydro_over_CH_significant"] is None
    f=field_coupling(1e-10,1e-9,1e-9,1e-20,1.,1e-10,1e-30,1.,POLICY)
    assert f["phi_coupling_current_relative"] is None
    assert f["mu_coupling_relative_diagnostic"]==pytest.approx(1e-4)
    with pytest.raises(ValueError): field_coupling(np.nan,0,0,0,1,0,0,1,POLICY)


def test_schedule_identity_cost_and_session():
    from sloshing.multiphase.benchmarks.phase_rate_series import schedules
    schedule=schedules()[0]
    assert schedule.sha256==L0_SHA and schedule.nsteps==1068
    rows=[dict(timing=dict(total_s=2.)) for _ in range(50)]
    r=cost_forecast(rows,1018,100.,120.)
    assert r["authorized"] and r["conservative_s_per_step"]==2.5
    assert not cost_forecast(rows,10180,100.,120.)["authorized"]
    assert session_result(None)=="complete" and session_result("error")=="failed"
