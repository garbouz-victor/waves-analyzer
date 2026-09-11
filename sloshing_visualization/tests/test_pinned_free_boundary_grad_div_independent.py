"""Independent grad-div algebra, tiny native audits and intentional corruptions.

Synthetic fields are algebra fixtures only, never PW2 trajectory evidence.
"""
from dataclasses import replace
import importlib.metadata
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest
from skfem import MeshTri, MeshTri2

from sloshing.pinned_wetting import free_boundary as production
from sloshing.pinned_wetting import free_boundary_run as execution
from sloshing.pinned_wetting import free_boundary_verify as audit
from sloshing.pinned_wetting import free_boundary_compare as compare

ROOT = Path(__file__).resolve().parents[2]


def _algebra_mesh(curved=False):
    linear = MeshTri.init_tensor(np.linspace(-1., 1., 4), np.linspace(-10., 0., 4))
    linear = linear.with_boundaries({"left": lambda x: x[0] == -1., "right": lambda x: x[0] == 1.,
        "bottom": lambda x: x[1] == -10., "surface": lambda x: x[1] == 0.})
    mesh = replace(MeshTri2.from_mesh(linear), _boundaries=linear.boundaries)
    if curved:
        points = mesh.p.copy()
        points[1] += .04*np.sin(np.pi*(points[0]+1.)/2.)*((points[1]+10.)/10.)**2
        mesh = replace(mesh, doflocs=points)
    return mesh


@pytest.mark.parametrize("curved", [False, True])
@pytest.mark.parametrize("gamma", [0., .1])
def test_independent_grad_div_all_rows_power_symmetry_and_PSD(curved, gamma):
    mesh = _algebra_mesh(curved)
    op = audit.rebuild(mesh.p, mesh.t, mesh.boundaries, 8, grad_div_gamma=gamma)
    rng = np.random.default_rng(819)
    velocity = rng.normal(scale=.2, size=op["velocity"].N)
    acceleration = rng.normal(scale=.3, size=len(velocity))
    pressure = rng.normal(scale=.4, size=op["pressure"].N)
    without = audit._actions(op, velocity, acceleration, pressure, 8)
    with_g = audit._actions(op, velocity, acceleration, pressure, 8, grad_div_gamma=gamma)
    G = op["grad_div"]
    np.testing.assert_allclose(with_g["momentum"]-without["momentum"], G@velocity, rtol=1e-12, atol=2e-14)
    np.testing.assert_allclose(with_g["weak_divergence"], without["weak_divergence"], rtol=0., atol=0.)
    np.testing.assert_allclose(with_g["grad_div_power"], velocity@(G@velocity), rtol=1e-12, atol=2e-14)
    np.testing.assert_allclose(G.toarray(), G.toarray().T, rtol=0., atol=2e-14)
    assert np.linalg.eigvalsh(G.toarray()).min() > -2e-12
    assert with_g["grad_div_power"] >= 0.
    assert with_g["bulk_power"] == without["bulk_power"]
    assert with_g["left_power"] == without["left_power"]
    # This is a separate production-vs-independent operator comparison; the
    # native verifier itself does not import the production operator.
    controls = production.Controls(nx=3, nz=3, grad_div_gamma_m2_s=gamma)
    implemented = production.Operators(mesh, controls).Kd
    np.testing.assert_allclose(implemented.toarray(), G.toarray(), rtol=2e-13, atol=2e-14)


def test_physical_and_augmented_free_tractions_cannot_be_conflated():
    mesh = _algebra_mesh()
    op = audit._basis(mesh.p, mesh.t, mesh.boundaries, 8)
    gamma = .1
    # v=.25 X has div v=.5, 2nu D=.005 I. p=.055 means physical
    # traction=-.05 n, whereas the AUGMENTED numerical traction is zero.
    velocity = np.zeros(op["velocity"].N)
    for axis, component in enumerate(op["components"]):
        velocity[component] = .25*mesh.p[axis]
    pressure = np.full(op["pressure"].N, .055)
    norms = audit._traction_norms(op, velocity, pressure, "surface", 8, False, gamma)
    np.testing.assert_allclose(norms["physical"], .05*np.sqrt(2.), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(norms["grad_div"], .05*np.sqrt(2.), rtol=1e-12, atol=1e-14)
    assert norms["augmented"] < 1e-13
    np.testing.assert_allclose(norms["grad_div_sampled_max"], .05, rtol=1e-12, atol=1e-14)
    left = audit._traction_norms(op, velocity, pressure, "left", 8, False, gamma)
    assert left["grad_div"] < 1e-14
    assert abs(left["augmented"]-left["physical"]) < 1e-14


@pytest.fixture(scope="module")
def grad_div_cases(tmp_path_factory):
    directory = tmp_path_factory.mktemp("independent_grad_div_tiny_native")
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "scikit-fem", "h5py")}
    cases = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(execution, "REL_RUNS", str(directory))
        patch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": "grad_div_unit_fixture"})
        for gamma in (0., .1):
            controls = production.Controls(nx=4, nz=4, dt=.0025, t_end=.01,
                hydrostatic_split=True, skew_divergence=True, rezone_interval_s=.005,
                rezone_strategy="straight_interior", grad_div_gamma_m2_s=gamma)
            result = execution.run_case(ROOT, controls, "grad-div-unit", {}, lambda *args: None)
            cases[gamma] = Path(result["native"])
    return cases


def _copy(native, tmp_path):
    destination = tmp_path/"mutated"
    shutil.copytree(native.parent, destination)
    return destination/"native.h5"


def _matrix(group, matrix):
    matrix = matrix.tocsr()
    for key in ("data", "indices", "indptr"):
        if key in group:
            del group[key]
        group.create_dataset(key, data=getattr(matrix, key))
    group.attrs["shape"] = matrix.shape


def _gamma_identity(path, gamma):
    with h5py.File(path, "r+") as h:
        identity = json.loads(h.attrs["numerical_identity"])
        identity["controls"]["grad_div_gamma_m2_s"] = gamma
        h.attrs["numerical_identity"] = json.dumps(identity)
        h.attrs["numerical_digest"] = audit._canonical(identity)
    (path.parent/"identity.json").write_text(json.dumps(identity))


def test_real_grad_div_native_reconstructs_numerical_work_but_keeps_raw_physical_gate(grad_div_cases):
    native = grad_div_cases[.1]
    report = audit.verify_native(native)
    assert report["passed"], report
    assert report["grad_div_gamma_m2_s"] == .1
    assert report["split_energy_per_density"]["signed_numerical_grad_div_work"] < 0.
    assert report["maxima"]["surface_grad_div_traction_L2_per_density"] > 0.
    assert report["maxima"]["saved_grad_div_work_discrepancy"] < 1e-12
    assert report["maxima"]["saved_cumulative_grad_div_work_discrepancy"] < 1e-12
    assert report["limits"]["continuous_energy_relative"] == .05
    with h5py.File(native, "r") as h:
        rows = [json.loads(h["states"][key].attrs["diagnostics"]) for key in sorted(h["states"])]
    raw = max(abs((row["energy"]-audit.E0+row["cumulative_bulk_loss"]+row["cumulative_left_wall_loss"])/audit.E0)
              for row in rows)
    np.testing.assert_allclose(report["maxima"]["continuous_energy_relative"], raw, rtol=0., atol=1e-11)
    assert audit.verify_native(grad_div_cases[0.])["passed"]


@pytest.mark.parametrize("mutation,expected", [
    ("missing_G", "missing_initial_grad_div_operator"),
    ("wrong_G", "initial_grad_div_operator_mismatch"),
    ("negative_G", "initial_grad_div_operator_mismatch"),
    ("asymmetric_G", "initial_grad_div_operator_not_symmetric"),
    ("wrong_work_sign", "numerical_grad_div_work_wrong_sign"),
    ("wrong_work_value", "limit:saved_grad_div_work_discrepancy"),
    ("wrong_cumulative_work", "limit:saved_cumulative_grad_div_work_discrepancy"),
    ("missing_work", "missing_numerical_grad_div_work"),
    ("wrong_gamma", "initial_grad_div_operator_mismatch"),
    ("negative_gamma", "invalid_numerical_grad_div_coefficient"),
])
def test_grad_div_native_mutations_are_rejected(grad_div_cases, tmp_path, mutation, expected):
    path = _copy(grad_div_cases[.1], tmp_path)
    if mutation in ("wrong_gamma", "negative_gamma"):
        _gamma_identity(path, .2 if mutation == "wrong_gamma" else -.1)
    else:
        with h5py.File(path, "r+") as h:
            if mutation == "missing_G":
                del h["initial_grad_div_operator"]
            elif mutation in ("wrong_G", "negative_G", "asymmetric_G"):
                G = audit._read_matrix(h["initial_grad_div_operator"])
                if mutation == "asymmetric_G":
                    G = G.tolil(); G[0, 1] += .01; G = G.tocsr()
                else:
                    G = G*(.5 if mutation == "wrong_G" else -1.)
                _matrix(h["initial_grad_div_operator"], G)
            else:
                group = h["states/00000002"]
                row = json.loads(group.attrs["diagnostics"])
                if mutation == "wrong_work_sign":
                    row["numerical_grad_div_work"] = 1e-5
                elif mutation == "wrong_work_value":
                    row["numerical_grad_div_work"] -= 1e-5
                elif mutation == "wrong_cumulative_work":
                    row["cumulative_numerical_grad_div_work"] -= 1e-5
                else:
                    del row["numerical_grad_div_work"]
                group.attrs["diagnostics"] = json.dumps(row)
    report = audit.verify_native(path)
    assert not report["passed"] and expected in report["failures"], report


def test_fields_computed_without_grad_div_cannot_be_relabelled_as_stabilized(grad_div_cases, tmp_path):
    path = _copy(grad_div_cases[0.], tmp_path)
    _gamma_identity(path, .1)
    with h5py.File(grad_div_cases[.1], "r") as correct, h5py.File(path, "r+") as wrong:
        _matrix(wrong["initial_grad_div_operator"], audit._read_matrix(correct["initial_grad_div_operator"]))
    report = audit.verify_native(path)
    assert not report["passed"] and "limit:momentum_residual" in report["failures"], report


def test_refinement_cannot_mix_gamma_and_publishes_numerical_loss_separately(grad_div_cases):
    with pytest.raises(ValueError, match="different_numerical_formulations"):
        compare.compare_runs(grad_div_cases[0.], grad_div_cases[.1], require_target=False)
    same = compare.compare_runs(grad_div_cases[.1], grad_div_cases[.1], require_target=False)
    assert same["passed"]
    differences = same["energy_observable_differences_per_density"]
    assert differences["numerical_grad_div_dissipation_rate"] == 0.
    assert same["accepted_strong_divergence_L2_difference"] == 0.
