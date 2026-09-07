"""Read-only full-operator check of archived converged Ritz vectors.

The initial run used a terminal Arnoldi recurrence estimate. This audit includes
actual floating-point operator/basis defects; previous output is not rewritten.
No time evolution, nonlinear solve or physical-parameter change occurs.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig
from sloshing.multiphase.benchmarks.spectral_ch import load_operator, RESULTS, read_json
from sloshing.multiphase.provenance import source_hash


def main():
    target = RESULTS/"spectrum/full_ritz_residual_audit.json"
    if target.exists():
        raise ValueError("Ritz audit output already exists")
    op, _, record = load_operator()
    spectrum = read_json(RESULTS/"spectrum/spectrum.json")
    with np.load(RESULTS/"krylov/basis.npz") as data:
        V = data["V"][:, :spectrum["dimension"]].copy()
        h = data["hessenberg"][:spectrum["dimension"]].copy()
    values, vectors = eig(h)
    rows = []
    for old in spectrum["ritz"]:
        if not (old["significant"] and old["converged"]):
            continue
        value = complex(old["real"], old["imag"])
        j = np.argmin(abs(values-value)); v = V @ vectors[:, j]
        residual = op.norm(op.apply(v)-value*v)/(max(abs(value), 1.)*op.norm(v))
        rows.append({"mode": old["mode"], "real": old["real"], "imag": old["imag"],
            "archived_recurrence_estimate": old["relative_Ritz_residual"],
            "full_operator_relative_residual": residual, "passed": residual <= 1e-6})
    result = {"status": "complete", "qualification_status": "passed" if rows and all(r["passed"] for r in rows) else "failed",
        "scope": "independent actual operator residual for all archived converged significant Ritz vectors",
        "count": len(rows), "max_full_operator_relative_residual": max(r["full_operator_relative_residual"] for r in rows),
        "rows": rows, "operator_fingerprint": op.fingerprint,
        "provenance": {**record["provenance"], "multiphase_source_sha256": source_hash(),
            "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("rows", "provenance")}))
    if result["qualification_status"] != "passed":
        raise SystemExit("STOP: full-operator Ritz residual failed")


if __name__ == "__main__":
    main()
