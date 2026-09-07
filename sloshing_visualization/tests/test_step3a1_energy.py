"""A decreasing E is not an energy budget. No FEM dependencies."""
from copy import deepcopy
from pathlib import Path
import numpy as np
import pytest
from sloshing.multiphase.energy_validation import validate_energy, read_history


def energy_history(multiplier=1.):
    rows=[]
    for t in np.linspace(0.,1.,11):
        rows.append({"time":t,"E_total":1.-.1*t,"E_interface":1.-.1*t,
            "E_kin":0.,"E_gravity":0.,"E_wall":0.,"mass_error_relative_to_domain":0.,
            "viscous_dissipation":0.,"CH_dissipation":.1*multiplier,"slip_dissipation":0.,
            "cumulative_viscous_dissipation":0.,"cumulative_CH_dissipation":.1*multiplier*t,
            "cumulative_slip_dissipation":0.})
    return rows


def test_monotone_energy_bogus_dissipation_cannot_pass():
    good=validate_energy(energy_history())
    assert good["qualified"] and good["absolute_budget_defect"]<1e-15
    bad=validate_energy(energy_history(1000.))
    assert bad["energy_growth_ok"]
    assert not bad["energy_budget_ok"] and not bad["qualified"]
    assert bad["relative_budget_defect_to_energy_change"]==pytest.approx(999.)


@pytest.mark.parametrize("key",["E_total","E_kin","E_gravity","E_interface","E_wall",
    "CH_dissipation","viscous_dissipation","slip_dissipation","cumulative_CH_dissipation",
    "energy_budget_relative"])
@pytest.mark.parametrize("value",[float("nan"),float("inf")])
def test_nonfinite_energy_fails_loudly(key,value):
    rows=energy_history();rows[2][key]=value
    with pytest.raises(ValueError,match="Nonfinite"):
        validate_energy(rows)


def test_saved_integral_is_checked_against_powers_not_trusted():
    rows=energy_history()
    rows[-1]["cumulative_CH_dissipation"]=0.
    assert not validate_energy(rows)["qualified"]


def test_component_scale_and_energy_change_floor_are_explicit():
    rows=energy_history()
    for row in rows:
        row["E_gravity"]=-10.;row["E_wall"]=10.
    result=validate_energy(rows)
    assert result["E_scale_initial_J_per_m"]==21.
    assert result["legacy_relative_budget_defect"]>=result["relative_budget_defect_to_initial_scale"]
    equilibrium=deepcopy(rows)
    for row in equilibrium:
        row["E_interface"]=row["E_total"]=1.
        row["CH_dissipation"]=row["cumulative_CH_dissipation"]=0.
    assert validate_energy(equilibrium)["relative_budget_defect_to_energy_change"] is None


def test_historical_456_J_run_is_rejected_without_changing_data():
    root=Path(__file__).resolve().parents[1]
    path=root/"validation_results/step3/contact/theta60_resolved/history.csv"
    result=validate_energy(read_history(path))
    assert result["energy_growth_ok"] and not result["energy_budget_ok"]
    ch=result["dissipation_components"]["CH"]
    assert ch["first_interval_J_per_m"]==pytest.approx(456.1060141331811)
    assert ch["largest_interval_fraction"]>.9999
    assert result["max_local_budget_defect_time_s"]==.01
    assert not result["qualified"]


def test_cancelled_initial_total_energy_does_not_corrupt_characteristic_gate():
    rows=energy_history()
    for row in rows:
        row["E_gravity"]=-1.
        row["E_total"]-=1.
    result=validate_energy(rows)
    assert result["qualified"]
    assert result["E_scale_initial_J_per_m"]==2.
    assert not result["legacy_metric_gates"]
    assert result["relative_budget_defect_to_initial_scale"]<1e-15
