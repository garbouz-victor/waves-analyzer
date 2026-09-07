from pathlib import Path
import numpy as np
import pytest
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.time_integrator import IntegrationSchedule,derivative_coefficients


def config(**changes):
    path=Path(__file__).resolve().parents[1]/"configs/step3/benchmark_flat.json"
    return ModelConfig.from_json(path).changed(**changes)


def integrate(c):
    schedule=IntegrationSchedule(c)
    values=[1.];times=[0.];phases=[];old=older=1.;spacing=None
    for step in range(schedule.nsteps):
        spec=schedule.step(step)
        if spec["scheme"]=="bdf2":
            assert spacing==spec["dt"]
        a,b,d=derivative_coefficients(2 if spec["scheme"]=="bdf2" else 1,spec["dt"],spec["scheme"])
        current=(-b*old-d*older)/(a+1.)  # y'=-y, no duplicated integrator formula
        older,old=old,current;spacing=spec["dt"]
        values.append(current);times.append(spec["time"]);phases.append(spec["phase"])
    return np.array(times),np.array(values),phases


def test_be_microsteps_then_correct_main_history():
    c=config(startup_dt=.001,startup_t_end=.01,dt=.02,t_end=.11)
    times,values,phases=integrate(c)
    assert phases[:10]==["startup_be"]*10
    assert phases[10]=="main_be_history" and phases[11]=="bdf2"
    assert times[9:13]==pytest.approx([.009,.01,.03,.05])
    assert values[10]==pytest.approx((1+.001)**-10)
    assert values[11]==pytest.approx(values[10]/1.02)


def test_bdf2_convergence_after_declared_restart():
    errors=[]
    for dt in (.04,.02,.01):
        c=config(dt=dt,t_end=.84,startup_dt=dt**2/4,startup_t_end=.04)
        times,values,_=integrate(c)
        errors.append(abs(values[-1]-np.exp(-times[-1])))
    assert errors[0]/errors[1]>3. and errors[1]/errors[2]>3.


def test_explicit_nonuniform_be_blocks_never_use_variable_bdf2():
    c=config(startup_stages=({"dt":.001,"t_end":.01},{"dt":.002,"t_end":.02}),dt=.01,t_end=.05)
    _,_,phases=integrate(c)
    assert phases[:15]==["startup_be"]*15 and phases[15]=="main_be_history"


def test_schedule_rejects_nonintegral_block():
    with pytest.raises(ValueError,match="integral"):
        IntegrationSchedule(config(startup_dt=.003,startup_t_end=.01))
