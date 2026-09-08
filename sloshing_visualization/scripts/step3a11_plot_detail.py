"""Post-measurement visual zoom only; no timing, classification, or policy change."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from step3a11_io_attribution import ROOT, read, sha, output, policy_guard

policy_guard()
source = ROOT/"benchmarks/bind_burst.json"
rows = read(source)["records"]
target = ROOT/"analysis/bind_append_phase_breakdown_ms.pdf"
if target.exists(): raise FileExistsError(target)
fig, ax = plt.subplots(figsize=(10, 4))
x = [r["repetition"] for r in rows]; bottom = np.zeros(len(rows))
for key, label in (("encode_hash_s", "prepare/hash/encode"), ("open_s", "open"),
                   ("write_s", "write"), ("flush_s", "flush"), ("fsync_s", "fsync"), ("close_s", "close")):
    values = np.array([r[key]*1000 for r in rows])
    ax.bar(x, values, bottom=bottom, width=.95, label=label); bottom += values
total = np.array([r["total_append_s"]*1000 for r in rows])
ax.bar(x, total-bottom, bottom=bottom, width=.95, color="lightgrey", label="probes/other overhead")
ax.set(xlabel="bind burst repetition", ylabel="append phase wall [ms]",
       title="Visual zoom only: 0/120 events reach the unchanged 1000 ms stall threshold")
ax.legend(ncol=4, fontsize=8); fig.tight_layout(); fig.savefig(target); plt.close(fig)
output(ROOT/"analysis/plot_detail_provenance.json", {
    "role": "post-measurement presentation only; all original samples and rules retained",
    "script_sha256": sha(__file__), "input_sha256": sha(source), "pdf_sha256": sha(target),
    "policy_sha256": sha(ROOT/"policy.json"), "original_full_scale_figure_unchanged": True,
    "new_benchmark_measurements": False})
