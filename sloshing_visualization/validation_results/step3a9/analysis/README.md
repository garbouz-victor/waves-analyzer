# Read-only prefix analysis — NOT full transfer qualification

The complete requested horizon was **not authorized by the frozen cost gate**.
Full-CHNS data here stop at accepted step127, time6.148107277665776e-10s.
`full_coupling_L0/COMPLETE.json` does not exist. No curve extrapolates full CHNS
to t=1e-4. In figure legends, “full L0” identifies the formulation/schedule,
not completion. `solver_iterations.pdf` also shows the prescribed later block
boundaries as reference lines; they are not executed intervals.

`prefix_metrics.json` uses the unchanged `scalar_coupling` / `field_coupling`
helpers, consistent historical P2 mass matrix, verified coefficient archives,
identical clocks, and only already-reached predetermined common checkpoints.
There is no interpolation, new threshold, Newton solve or state mutation.
The final-time transfer gate function is deliberately not called on a prefix.

The plots are rendered with the existing frozen plotting functions. For the
coupling plots both histories are restricted to the same measured prefix.
`energy_budget.pdf` is a byte-identical filename alias of `budget.pdf`.
The `strong_divergence_L0.pdf` and `weak_continuity_L0.pdf` filenames identify
the incomplete L0 dataset; their statistics and curves are prefix-only.

Post-run editorial additions to summary.json/csv and STEP3A9_REPORT.md spell
out the cost STOP and final local tests. They do not alter any frozen source,
gate, measured history, checkpoint or provenance. Re-running the generic
report summarizer regenerates its shorter, still-unqualified verdict.
