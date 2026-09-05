"""Small reproducible CSV/JSON + scientific PDF/PNG output, no animation."""

import csv
import json
import os
from pathlib import Path
import tempfile


def write_tables(output, name, rows, summary=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (output/f"{name}.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    with (output/f"{name}.json").open("w", encoding="utf-8") as handle:
        json.dump({"rows": rows, "summary": summary or {}}, handle, indent=2, allow_nan=False)


def pyplot():
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="sloshing-mpl-"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_figure(figure, output, name):
    for extension in ("pdf", "png"):
        figure.savefig(Path(output)/f"{name}.{extension}", dpi=180, bbox_inches="tight")
