# Independent coordinate-portability addendum

This addendum does not replace or modify the original signed `review.md`,
`independent_review.json`, `native_audit.json`, their exact auditor snapshot,
or any native/production numerical source. It reviews a verifier-only correction.

## Measured defect and narrow correction

The original native trajectories were generated with NumPy 1.21.5 / SciPy 1.7.3.
Reconstructing the same tanh-graded mesh with NumPy 2.0.2 / SciPy 1.13.1 changed
surface x coordinates by at most 2.220446049250313e-16 m, bulk x coordinates by
at most 1.1102230246251565e-16 m, and bulk z coordinates by exactly zero.
These values were measured independently for all three preserved native files.
Shapes and triangle connectivity matched exactly. Reconstruction with the
original runtime reproduced every coordinate exactly.

The previous coordinate `array_equal` therefore rejected elementary-function
and P2 midpoint roundoff, not a physically different mesh. Only coordinate
comparison now permits an absolute bound of 8 binary64 machine epsilons times
the physical extent of each axis: 1.7763568394002505e-15 m horizontally and
1.7763568394002505e-14 m vertically. This accommodates a few arithmetic and
elementary-function rounding errors; it is not a mesh-convergence tolerance.
There is no default relative `allclose` tolerance. Shape and finiteness checks,
triangle connectivity, free/fixed DOFs and operator support remain strict.
Physical, residual, endpoint, coating and exact-slope gates are unchanged.

## Actual clean-runtime verification

One kernel-locked managed job, session43564, completed with exit0. In the fresh
dependency-only environment the independent verifier suite passed 21 tests,
zero failures/errors/skips, pytest elapsed8.21s. Added checks accept one-ULP
coordinate differences but reject a real native surface/bulk displacement of
1e-8m, modified triangle connectivity, shape changes and nonfinite coordinates.

The same managed job then re-audited every SDIRK2 stage and operator of the
unchanged actual baseline pilot in the clean runtime. Equation/BC/source/coating/
energy checks passed; the sole native failure remained `limit:max_slope`.
Measured job-body time was14.0716s. This is a successful portability correction,
not a successful physical qualification or a substitute for the forthcoming
fresh three-native package audit. The scientific blocker and target NOT_RUN
verdict remain unchanged.

New verifier SHA256:
`c173d7e0c6354de10e88aa6718395fdf13a50dc120add567072b2ab072c211c0`.
Its exact snapshot is `native_auditor_portability_source.py` in corrected output.

Evidence in `sloshing_visualization/output/pinned_wetting/left_navier_right_noslip`:

- `mesh_portability_numpy_2.0.2.json`, SHA256
  `6caf2489124f370d9b8d898fe0f3cf84121eabe273e30561591184431c0f3f13`;
- `mesh_portability_numpy_1.21.5.json`, SHA256
  `a9be9e52a957ac151fa4c51f13f422fa0d71536039bf04ba3cd1142599e6dc40`;
- `portability_native_probe.json`, SHA256
  `e33fe5c7aa95b49445a9e7b14caa2fb2858da8542d6c9ff2b9dd58dd7568ea0c`.

Test evidence is adjacent to this addendum: `reviewer_portability_tests.xml`,
SHA256 `9646dd17e1a422f2a4b0d96d10e599dfadd2f15c4eef94a1401f61257a49fd39`.
Reviewer context: `/root/physics_review`. No old signed evidence was rewritten.
