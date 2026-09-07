"""Re-evaluate frozen nonlinear initial state only: NO time advance."""
import json
import numpy as np
from dolfinx import fem
from sloshing.multiphase.benchmarks.spectral_ch import scientific_input, reserve, write_json
from sloshing.multiphase.linearized_ch import array_hash
from sloshing.multiphase.diagnostics import Diagnostics, interface_resolution


def main():
    output = reserve("linear_operator/initial_audit")
    s, eq, psi, coefficients, provenance = scientific_input()
    equilibrium = Diagnostics(s).measure()
    initial = fem.Function(eq.function_space)
    initial.x.array[:] = coefficients; initial.x.scatter_forward()
    s.initialize(initial)  # Existing consistent weak chemical-potential initialization.
    row = Diagnostics(s).measure()
    resolution = interface_resolution(s)
    actual_phi = s.state.sub(2).collapse().x.array
    if not np.array_equal(coefficients, actual_phi) or s.step_number != 0 or s.time != 0:
        raise RuntimeError("Initial field bytes changed or unexpected time advance")
    result = {"status": "complete", "scope": "nonlinear initial-state audit; zero accepted intervals",
        "provenance": provenance, "initial": row, "equilibrium": equilibrium,
        "initial_excess_energy": row["E_total"]-equilibrium["E_total"],
        "tau_E": (row["E_total"]-equilibrium["E_total"])/row["CH_dissipation"],
        "resolution": resolution, "actual_initial_phi_sha256": array_hash(actual_phi),
        "psi_sha256": array_hash(psi), "accepted_intervals": 0,
        "nonlinear_transient_qualification": "not_run; cost STOP remains in force"}
    write_json(output/"status.json", result)
    print(json.dumps({"D_CH_0": row["CH_dissipation"], "excess_energy": result["initial_excess_energy"],
        "tau_E": result["tau_E"], "actual_initial_phi_sha256": result["actual_initial_phi_sha256"],
        "accepted_intervals": 0}))


if __name__ == "__main__":
    main()
