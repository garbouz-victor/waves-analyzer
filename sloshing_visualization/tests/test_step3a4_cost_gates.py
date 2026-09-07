import pytest
from sloshing.multiphase.ch_validation_policy import (AnalysisBudget, ValidationPolicy,
    run_sparse_analysis, authorizations, linear_gates, LINEAR_LIMITS)
from sloshing.multiphase.performance import Performance


def test_expensive_full_CHNS_does_not_skip_affordable_sparse_analysis():
    cost = 1796*9.3186596433
    assert cost > 3*3600 and 7*cost > 10*3600
    executed = []
    result = run_sparse_analysis(lambda: executed.append(True) or {"passed": True},
        reduced_pass=True, analysis_budget=AnalysisBudget())
    assert executed == [True] and result["passed"]


def test_analysis_has_its_own_finite_budget():
    with pytest.raises(RuntimeError, match="analysis cost gate"):
        AnalysisBudget(7201.).check()
    with pytest.raises(RuntimeError, match="analysis cost gate"):
        AnalysisBudget().check(3601.)
    with pytest.raises(ValueError, match="Reduced"):
        run_sparse_analysis(lambda: pytest.fail("must not execute"), reduced_pass=False,
                            analysis_budget=AnalysisBudget())


def test_separate_scientific_roles_and_no_threshold_relaxation():
    flags = authorizations(reduced=True, sparse=True)
    assert flags["full_sparse_linear_verified"]
    assert not flags["isolated_CH_execution_authorized"]
    assert not flags["full_CHNS_probe_authorized"]
    assert "nonlinear_execution_authorized" not in flags
    assert flags["full_model_temporal_qualification"] == "not_qualified"
    assert all(linear_gates(LINEAR_LIMITS).values())
    assert not all(linear_gates({**LINEAR_LIMITS, "integrated_D2_relative_error": .01001}).values())
    assert not all(linear_gates(LINEAR_LIMITS, full_sparse=True).values())
    assert ValidationPolicy().certified_cells == 8


def test_performance_reports_unmeasured_as_null_and_measures_errors():
    perf = Performance()
    with pytest.raises(ValueError), perf.measure("linear_solve", dt=.1):
        raise ValueError("measurement persists on failure")
    value = perf.snapshot()
    assert value["categories"]["linear_solve"]["count"] == 1
    assert value["categories"]["nonlinear_SNES"] is None
    assert value["records"]["linear_solve"][0]["dt"] == .1
    assert value["peak_RSS_bytes"] > 0


def test_step3a4_cannot_write_historical_paths(tmp_path, monkeypatch):
    from sloshing.multiphase.benchmarks import cost_efficient_ch as stages
    monkeypatch.setattr(stages, "RESULTS", tmp_path/"step3a4")
    for path in ("../step3a3", str(tmp_path/"outside"), "."):
        with pytest.raises(ValueError, match="Historical/external"):
            stages.reserve(path)


def test_legacy_planner_entrypoint_also_runs_sparse_when_full_CHNS_too_expensive(tmp_path, monkeypatch):
    import numpy as np
    from sloshing.multiphase.benchmarks import spectral_ch
    from sloshing.multiphase import ch_timestep_planner as planner
    from sloshing.multiphase.linearized_ch import LinearizedCH
    op = LinearizedCH(np.eye(3), np.diag([1., 100., 10000.]), np.eye(3))
    metadata = {"provenance": {"equilibrium_fingerprint": "eq", "perturbation_coefficient_sha256": "psi"},
        "CH_dominated_justification": {"cost_seconds_per_step_conservative": 9.3186596433}}
    reference = {"Krylov_sha256": "krylov", "rho_relevant": 10000., "dimension": 3}
    basis = {"V": np.eye(3), "projected_L": -np.diag([1., 100., 10000.]), "beta": 1.}
    monkeypatch.setattr(spectral_ch, "qualified_reference", lambda: (op, None, metadata, reference, basis, None))
    monkeypatch.setattr(spectral_ch, "read_json", lambda path: {"qualification_status": "passed", "mode": {
        "operator_fingerprint": op.fingerprint, "Krylov_sha256": "krylov"}})
    monkeypatch.setattr(spectral_ch, "reserve", lambda name: tmp_path)
    monkeypatch.setattr(planner, "make_plan", lambda *a, **kw: {"qualification_status": "linear_reduced_pass",
        "total_steps": 1796, "blocks": [], "predicted_errors": {k: 0. for k in LINEAR_LIMITS}, "dt0": 1e-12})
    calls = []
    def verify(*args, **kwargs):
        calls.append(True)
        return {**{k: 0. for k in LINEAR_LIMITS}, "max_mass_domain": 0.}
    monkeypatch.setattr(planner, "verify_full_sparse_schedule", verify)
    result = spectral_ch.planner_stage()
    assert calls == [True]
    assert result["gates"]["full_sparse_BE"]
    assert not result["gates"]["cost"]
