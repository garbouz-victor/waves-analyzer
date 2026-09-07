import numpy as np
import pytest
from sloshing.multiphase.settling import validate_settling,physical_window_rate


def series(dt,rate=0.):
    times=np.linspace(0.,.2,round(.2/dt)+1)
    observations=[];history=[]
    for t in times:
        observations.append({"time":t,"fits":[{"theta_deg":60+rate*t}]*3,
            "crossings":[{"coordinate":-.4},{"coordinate":.4}]})
        history.append({"time":t,"E_kin":0.,"speed_max_dof_sample":0.,"relative_mu_variation":0.})
    return observations,history


def test_identical_physical_rate_independent_of_dt():
    results=[validate_settling(*series(dt,.2),1.) for dt in (.01,.001)]
    assert all(r["qualified"] for r in results)
    assert [r["angle_rate_deg_per_s"] for r in results]==pytest.approx([.2,.2])


def test_tiny_change_per_step_with_large_rate_fails():
    result=validate_settling(*series(.0001,5.),1.)
    assert 5*.0001<.05  # the old per-step gate would pass
    assert not result["qualified"]
    assert result["angle_rate_deg_per_s"]==pytest.approx(5.)


def test_settled_and_insufficient_window():
    assert validate_settling(*series(.01),1.)["qualified"]
    obs,hist=series(.01)
    assert not validate_settling(obs[:5],hist[:5],1.)["qualified"]


def test_time_weighted_rate_with_irregular_sampling():
    times=np.array([0.,.00001,.001,.003,.1,.2])
    assert physical_window_rate(times,2+3*times,.1)==pytest.approx(3.)


def test_stationary_angle_alone_does_not_pass():
    obs,hist=series(.01)
    for row in hist:
        row["speed_max_dof_sample"]=.01
    assert not validate_settling(obs,hist,1.)["qualified"]
