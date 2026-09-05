# STEP 1.5 — complete test results

Measured on 2026-09-06: **73 passed, 9 warnings in 107.55 s**.
All 43 original STEP 1 tests remain passing. Default selection: 71;
`validation` marker: 2 (spatial orders and temporal contamination).

Command (both selections together):

```bash
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m pytest -v -m 'validation or not validation' \
  --junitxml=validation_results/tests.xml
```

Machine-readable evidence: [tests.xml](tests.xml).

## By file

| File | Passed |
|---|---:|
| [tests/test_config.py](../tests/test_config.py) | 11 |
| [tests/test_contact_comparison.py](../tests/test_contact_comparison.py) | 1 |
| [tests/test_contact_points.py](../tests/test_contact_points.py) | 3 |
| [tests/test_discrete_mode.py](../tests/test_discrete_mode.py) | 1 |
| [tests/test_divergence.py](../tests/test_divergence.py) | 2 |
| [tests/test_energy_decay.py](../tests/test_energy_decay.py) | 10 |
| [tests/test_integrator_consistency.py](../tests/test_integrator_consistency.py) | 9 |
| [tests/test_linear_mode.py](../tests/test_linear_mode.py) | 4 |
| [tests/test_manufactured_spatial.py](../tests/test_manufactured_spatial.py) | 4 |
| [tests/test_manufactured_symbolic.py](../tests/test_manufactured_symbolic.py) | 4 |
| [tests/test_no_slip.py](../tests/test_no_slip.py) | 5 |
| [tests/test_sdirk2_fem_mode.py](../tests/test_sdirk2_fem_mode.py) | 1 |
| [tests/test_sdirk2_scalar.py](../tests/test_sdirk2_scalar.py) | 9 |
| [tests/test_storage.py](../tests/test_storage.py) | 3 |
| [tests/test_symmetry.py](../tests/test_symmetry.py) | 3 |
| [tests/test_validation_cache.py](../tests/test_validation_cache.py) | 1 |
| [tests/test_volume.py](../tests/test_volume.py) | 1 |
| [tests/test_vorticity_regions.py](../tests/test_vorticity_regions.py) | 1 |

## Every collected case

| Test | Result |
|---|---|
| `tests/test_config.py::test_reject_invalid_config[changes0]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes1]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes2]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes3]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes4]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes5]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes6]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes7]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes8]` | PASS |
| `tests/test_config.py::test_reject_invalid_config[changes9]` | PASS |
| `tests/test_config.py::test_config_roundtrip_and_stable_run_id` | PASS |
| `tests/test_contact_comparison.py::test_p2_evaluation_and_cross_mesh_quadrature` | PASS |
| `tests/test_contact_points.py::test_contact_points_are_exactly_pinned` | PASS |
| `tests/test_contact_points.py::test_exact_quadratic_surface_slope` | PASS |
| `tests/test_contact_points.py::test_warning_uses_dimensional_eta_slope_without_extra_amplitude_factor` | PASS |
| `tests/test_discrete_mode.py::test_midpoint_against_independent_viscous_semidiscrete_eigenmode` | PASS |
| `tests/test_divergence.py::test_weak_and_strong_divergence_are_reported_separately` | PASS |
| `tests/test_divergence.py::test_actual_divergence_decreases_on_refinement` | PASS |
| `tests/test_energy_decay.py::test_energy_budget_on_every_internal_step` | PASS |
| `tests/test_energy_decay.py::test_zero_angle_stays_at_rest` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[left_u_max]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[right_w_max]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[bottom_u_max]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[volume_error]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[weak_divergence_l2]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[contact_error]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[energy_balance_relative]` | PASS |
| `tests/test_energy_decay.py::test_invalid_diagnostics_fail_loudly[max_step_energy_increase]` | PASS |
| `tests/test_integrator_consistency.py::test_source_stage_times_and_forced_energy_identity[midpoint]` | PASS |
| `tests/test_integrator_consistency.py::test_source_stage_times_and_forced_energy_identity[sdirk2]` | PASS |
| `tests/test_integrator_consistency.py::test_sdirk_stages_satisfy_monolithic_dae` | PASS |
| `tests/test_integrator_consistency.py::test_sdirk_production_invariants_and_rk_budget[0.01]` | PASS |
| `tests/test_integrator_consistency.py::test_sdirk_production_invariants_and_rk_budget[0.1]` | PASS |
| `tests/test_integrator_consistency.py::test_sdirk_production_invariants_and_rk_budget[1.0]` | PASS |
| `tests/test_integrator_consistency.py::test_sdirk_manufactured_smoke` | PASS |
| `tests/test_integrator_consistency.py::test_mms_temporal_second_order_on_fixed_mesh[midpoint]` | PASS |
| `tests/test_integrator_consistency.py::test_mms_temporal_second_order_on_fixed_mesh[sdirk2]` | PASS |
| `tests/test_linear_mode.py::test_release_direction_is_computed` | PASS |
| `tests/test_linear_mode.py::test_linear_amplitude_and_sign_response` | PASS |
| `tests/test_linear_mode.py::test_constant_surface_fixes_pressure_gauge_without_spurious_motion` | PASS |
| `tests/test_linear_mode.py::test_second_order_time_convergence` | PASS |
| `tests/test_manufactured_spatial.py::test_manufactured_spatial_smoke_in_default_suite` | PASS |
| `tests/test_manufactured_spatial.py::test_production_cannot_enable_manufactured_forcing_via_config` | PASS |
| `tests/test_manufactured_spatial.py::test_observed_manufactured_spatial_orders` | PASS |
| `tests/test_manufactured_spatial.py::test_mms_spatial_errors_are_not_time_error` | PASS |
| `tests/test_manufactured_symbolic.py::test_symbolic_incompressibility_and_wall_no_slip` | PASS |
| `tests/test_manufactured_symbolic.py::test_symbolic_kinematics_and_pinned_endpoints` | PASS |
| `tests/test_manufactured_symbolic.py::test_symbolic_momentum_including_stress_divergence` | PASS |
| `tests/test_manufactured_symbolic.py::test_symbolic_both_traction_components` | PASS |
| `tests/test_no_slip.py::test_all_six_wall_residuals_and_independent_boundary_probes` | PASS |
| `tests/test_no_slip.py::test_every_viscosity_preset_obeys_constraints[0.001]` | PASS |
| `tests/test_no_slip.py::test_every_viscosity_preset_obeys_constraints[0.01]` | PASS |
| `tests/test_no_slip.py::test_every_viscosity_preset_obeys_constraints[0.1]` | PASS |
| `tests/test_no_slip.py::test_every_viscosity_preset_obeys_constraints[1.0]` | PASS |
| `tests/test_sdirk2_fem_mode.py::test_fast_constrained_fem_mode_distinguishes_integrators` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[midpoint--1.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[midpoint--100.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[midpoint--10000.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--1.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--100.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--10000.0]` | PASS |
| `tests/test_sdirk2_scalar.py::test_second_order_scalar_decay` | PASS |
| `tests/test_sdirk2_scalar.py::test_l_stability_is_not_a_stability` | PASS |
| `tests/test_sdirk2_scalar.py::test_a_stability_on_imaginary_boundary_and_infinity` | PASS |
| `tests/test_storage.py::test_vorticity_sign_and_exact_fem_derivative` | PASS |
| `tests/test_storage.py::test_hdf5_roundtrip_and_csv` | PASS |
| `tests/test_storage.py::test_failed_write_is_marked_and_closed` | PASS |
| `tests/test_symmetry.py::test_mesh_and_release_preserve_reflection_symmetry[8]` | PASS |
| `tests/test_symmetry.py::test_mesh_and_release_preserve_reflection_symmetry[9]` | PASS |
| `tests/test_symmetry.py::test_wall_max_includes_quadratic_extrema_between_dofs` | PASS |
| `tests/test_validation_cache.py::test_interrupted_cache_is_neither_reused_nor_overwritten` | PASS |
| `tests/test_volume.py::test_volume_is_preserved_including_endpoint_basis_functions` | PASS |
| `tests/test_vorticity_regions.py::test_exact_clipped_region_integrals` | PASS |

## Visible warnings

One warning: original contact test reaches max slope 0.3037.
Eight warnings: small test meshes under-resolve strong divergence;
relative norms 0.173, 0.227, 0.260 and 0.465. These are deliberately
reported, not suppressed or replaced by the near-zero weak residual.
The medium physical runs have separate measured diagnostics in
[STEP1_5_REPORT.md](../STEP1_5_REPORT.md).
