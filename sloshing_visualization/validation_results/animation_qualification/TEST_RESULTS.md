# STEP 1.6 final local test results

Measured on the final STEP 1.6 source, Python 3.9.12 / tested dependencies,
with OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1, OMP_NUM_THREADS=1.
GitHub Actions is configured for Python 3.11 but has not been executed remotely.

| Invocation | Result | Wall time |
|---|---|---:|
| `python -m pytest -q` | 94 passed, 2 deselected, 14 warnings | 25.06 s |
| `python -m pytest -q -m validation` | 2 passed, 94 deselected | 91.17 s |
| Explicit ten STEP 1.6 test modules below | 23 passed, 5 warnings | 24.13 s |

There are **96 distinct passing tests**, including the 73 existing tests and
23 new tests. The last invocation is a subset, not 23 additional test cases.
Warnings are the existing pinned-slope breakdown test and intentionally
under-resolved tiny/coarse test meshes; they are not suppressed. There are
no failed or skipped tests. Full qualification runs are not unit tests.

All medium/fine short/full phases, full exact-overlay field comparison,
particle qualification (including a separate dense-snapshot PDE run), plots
and final summary completed. Rebuilding the summary from the same complete
HDF5 files produced exact equality of parsed JSON values. All nine required
PDF/PNG pairs and README/report local links were checked. HDF5, full raw
curves, trajectory NPZ and test XML are ignored generated artifacts.

Scientific decision remains **NOT QUALIFIED FOR PHYSICAL ANIMATION** at
alpha=0.02 degrees because the measured fine contact slope exceeds 0.3.
Passing tests does not override that model-qualification failure.

## Complete test inventory

Every entry below passed; `(validation)` marks the two opt-in expensive tests.

- `tests.test_amplitude_scaling::test_real_production_path_scales_with_tan_alpha` — PASS
- `tests.test_animation_candidate_config::test_animation_candidates_and_historical_config` — PASS
- `tests.test_config::test_config_roundtrip_and_stable_run_id` — PASS
- `tests.test_config::test_reject_invalid_config[changes0]` — PASS
- `tests.test_config::test_reject_invalid_config[changes1]` — PASS
- `tests.test_config::test_reject_invalid_config[changes2]` — PASS
- `tests.test_config::test_reject_invalid_config[changes3]` — PASS
- `tests.test_config::test_reject_invalid_config[changes4]` — PASS
- `tests.test_config::test_reject_invalid_config[changes5]` — PASS
- `tests.test_config::test_reject_invalid_config[changes6]` — PASS
- `tests.test_config::test_reject_invalid_config[changes7]` — PASS
- `tests.test_config::test_reject_invalid_config[changes8]` — PASS
- `tests.test_config::test_reject_invalid_config[changes9]` — PASS
- `tests.test_contact_comparison::test_p2_evaluation_and_cross_mesh_quadrature` — PASS
- `tests.test_contact_points::test_contact_points_are_exactly_pinned` — PASS
- `tests.test_contact_points::test_exact_quadratic_surface_slope` — PASS
- `tests.test_contact_points::test_warning_uses_dimensional_eta_slope_without_extra_amplitude_factor` — PASS
- `tests.test_discrete_mode::test_midpoint_against_independent_viscous_semidiscrete_eigenmode` — PASS
- `tests.test_divergence::test_actual_divergence_decreases_on_refinement` — PASS
- `tests.test_divergence::test_weak_and_strong_divergence_are_reported_separately` — PASS
- `tests.test_energy_decay::test_energy_budget_on_every_internal_step` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[bottom_u_max]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[contact_error]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[energy_balance_relative]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[left_u_max]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[max_step_energy_increase]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[right_w_max]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[volume_error]` — PASS
- `tests.test_energy_decay::test_invalid_diagnostics_fail_loudly[weak_divergence_l2]` — PASS
- `tests.test_energy_decay::test_zero_angle_stays_at_rest` — PASS
- `tests.test_integrator_consistency::test_mms_temporal_second_order_on_fixed_mesh[midpoint]` — PASS
- `tests.test_integrator_consistency::test_mms_temporal_second_order_on_fixed_mesh[sdirk2]` — PASS
- `tests.test_integrator_consistency::test_sdirk_manufactured_smoke` — PASS
- `tests.test_integrator_consistency::test_sdirk_production_invariants_and_rk_budget[0.01]` — PASS
- `tests.test_integrator_consistency::test_sdirk_production_invariants_and_rk_budget[0.1]` — PASS
- `tests.test_integrator_consistency::test_sdirk_production_invariants_and_rk_budget[1.0]` — PASS
- `tests.test_integrator_consistency::test_sdirk_stages_satisfy_monolithic_dae` — PASS
- `tests.test_integrator_consistency::test_source_stage_times_and_forced_energy_identity[midpoint]` — PASS
- `tests.test_integrator_consistency::test_source_stage_times_and_forced_energy_identity[sdirk2]` — PASS
- `tests.test_linear_mode::test_constant_surface_fixes_pressure_gauge_without_spurious_motion` — PASS
- `tests.test_linear_mode::test_linear_amplitude_and_sign_response` — PASS
- `tests.test_linear_mode::test_release_direction_is_computed` — PASS
- `tests.test_linear_mode::test_second_order_time_convergence` — PASS
- `tests.test_manufactured_spatial::test_manufactured_spatial_smoke_in_default_suite` — PASS
- `tests.test_manufactured_spatial::test_mms_spatial_errors_are_not_time_error` (validation) — PASS
- `tests.test_manufactured_spatial::test_observed_manufactured_spatial_orders` (validation) — PASS
- `tests.test_manufactured_spatial::test_production_cannot_enable_manufactured_forcing_via_config` — PASS
- `tests.test_manufactured_symbolic::test_symbolic_both_traction_components` — PASS
- `tests.test_manufactured_symbolic::test_symbolic_incompressibility_and_wall_no_slip` — PASS
- `tests.test_manufactured_symbolic::test_symbolic_kinematics_and_pinned_endpoints` — PASS
- `tests.test_manufactured_symbolic::test_symbolic_momentum_including_stress_divergence` — PASS
- `tests.test_no_slip::test_all_six_wall_residuals_and_independent_boundary_probes` — PASS
- `tests.test_no_slip::test_every_viscosity_preset_obeys_constraints[0.001]` — PASS
- `tests.test_no_slip::test_every_viscosity_preset_obeys_constraints[0.01]` — PASS
- `tests.test_no_slip::test_every_viscosity_preset_obeys_constraints[0.1]` — PASS
- `tests.test_no_slip::test_every_viscosity_preset_obeys_constraints[1.0]` — PASS
- `tests.test_particle_interpolation::test_actual_p2_interpolation_matches_fem_probes[4]` — PASS
- `tests.test_particle_interpolation::test_actual_p2_interpolation_matches_fem_probes[5]` — PASS
- `tests.test_particle_interpolation::test_overlay_integrates_both_p2_fields_without_visualization_sampling` — PASS
- `tests.test_particle_interpolation::test_overlay_resolves_discontinuous_vorticity_across_both_meshes` — PASS
- `tests.test_particle_snapshot_refinement::test_invalid_surface_crossing_is_stopped_not_reflected` — PASS
- `tests.test_particle_snapshot_refinement::test_linear_snapshot_interpolation_converges_second_order` — PASS
- `tests.test_particle_snapshot_refinement::test_particle_rk4_has_fourth_order_for_exact_fem_linear_flow` — PASS
- `tests.test_particle_wall_velocity::test_fem_particle_evaluation_preserves_wall_zero_and_limit` — PASS
- `tests.test_qualification_cli::test_standalone_particles_cannot_bypass_full_field_gate` — PASS
- `tests.test_qualification_comparison::test_cross_mesh_comparison_and_compatible_reuse` — PASS
- `tests.test_qualification_dataset::test_complete_prefix_continuation_matches_uninterrupted` — PASS
- `tests.test_qualification_dataset::test_exact_regional_slopes_include_cut_inside_p2_edge` — PASS
- `tests.test_qualification_dataset::test_prefix_target_cannot_be_overwritten_by_concurrent_creation` — PASS
- `tests.test_qualification_dataset::test_regional_derivatives_against_analytic_polynomial_integral` — PASS
- `tests.test_qualification_dataset::test_static_surface_interpretation_matches_actual_trace_nullspace` — PASS
- `tests.test_qualification_dataset::test_stiffly_accurate_stage_pressure_is_instantaneous_pressure[0.01]` — PASS
- `tests.test_qualification_dataset::test_stiffly_accurate_stage_pressure_is_instantaneous_pressure[1.0]` — PASS
- `tests.test_qualification_policy::test_divergence_and_particles_can_fail_otherwise_good_dataset` — PASS
- `tests.test_qualification_policy::test_slope_policy_cannot_be_overruled_by_passing_invariants` — PASS
- `tests.test_qualification_signals::test_horizontal_fem_line_quadrature_is_exact_for_p2_squared` — PASS
- `tests.test_qualification_signals::test_modal_event_estimator_on_independent_damped_signal` — PASS
- `tests.test_sdirk2_fem_mode::test_fast_constrained_fem_mode_distinguishes_integrators` — PASS
- `tests.test_sdirk2_scalar::test_a_stability_on_imaginary_boundary_and_infinity` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[midpoint--1.0]` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[midpoint--100.0]` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[midpoint--10000.0]` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--1.0]` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--100.0]` — PASS
- `tests.test_sdirk2_scalar::test_actual_scalar_stages_match_analytic_stability_function[sdirk2--10000.0]` — PASS
- `tests.test_sdirk2_scalar::test_l_stability_is_not_a_stability` — PASS
- `tests.test_sdirk2_scalar::test_second_order_scalar_decay` — PASS
- `tests.test_storage::test_failed_write_is_marked_and_closed` — PASS
- `tests.test_storage::test_hdf5_roundtrip_and_csv` — PASS
- `tests.test_storage::test_vorticity_sign_and_exact_fem_derivative` — PASS
- `tests.test_symmetry::test_mesh_and_release_preserve_reflection_symmetry[8]` — PASS
- `tests.test_symmetry::test_mesh_and_release_preserve_reflection_symmetry[9]` — PASS
- `tests.test_symmetry::test_wall_max_includes_quadratic_extrema_between_dofs` — PASS
- `tests.test_validation_cache::test_interrupted_cache_is_neither_reused_nor_overwritten` — PASS
- `tests.test_volume::test_volume_is_preserved_including_endpoint_basis_functions` — PASS
- `tests.test_vorticity_regions::test_exact_clipped_region_integrals` — PASS
