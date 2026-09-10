"""Genuine tiny asymmetric trajectories and explicit structural corruptions."""
from dataclasses import asdict
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest
from scipy.sparse import coo_matrix
from skfem import asm

from sloshing.config import SimulationConfig
from sloshing.fem_spaces import scalar_mass
from sloshing.pinned_wetting import asymmetric_run
from sloshing.pinned_wetting.asymmetric import LeftNavierConfig, LeftNavierDiagnostics, LeftNavierFEM
from sloshing.pinned_wetting.asymmetric_verify import (
    CONTACT, LABEL, MODEL_ID, SOURCE_KIND, _hash, _operators, _read_wall_matrix,
    _coordinate_check, channel_reference, comparison_metrics, verify_native,
)
from sloshing.pinned_wetting.mission import write_json

ROOT=Path(__file__).resolve().parents[2]


def make_case(directory, *, slip=.5, t_end=.02):
    directory.mkdir()
    c=LeftNavierConfig(nx=8,nz=12,dt=.005,t_end=t_end,snapshot_dt=.01,integrator="sdirk2",left_slip_length_m=slip)
    source_hash,source_files=asymmetric_run.numerical_identity(ROOT)
    physical=asymmetric_run.physical_identity(c)
    identity={"model_id":MODEL_ID,"run_id":directory.name,"source_hash":source_hash,"source_files":source_files,
        "physical":physical,"physical_config_hash":_hash(physical),"physical_model_label":LABEL,
        "source_kind":SOURCE_KIND,"synthetic":False,"role":"unit_test_real_asymmetric_trajectory"}
    write_json(directory/"identity.json",identity)
    write_json(directory/"resolved_case.json",{"model_id":MODEL_ID,"controls":asdict(c),"physical":physical,
        "physical_model_label":LABEL,"film_thickness_m":None,"retained_layer_volume_feedback":"neglected_at_leading_order",
        "extension_contract":{"record_height_tolerance_m":.0001}})
    for name in source_files:
        target=directory/"source_snapshot"/name
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/"sloshing_visualization/src/sloshing"/name,target)
    asymmetric_run.execute_case(ROOT,directory,c,identity,{},lambda *args:None)
    return directory/"native.h5"


def replace_wall(h,matrix):
    group=h["operators/K_wall_full"]
    matrix=matrix.tocsr()
    matrix.sum_duplicates(); matrix.eliminate_zeros(); matrix.sort_indices()
    for name in ("data","indices","indptr"):
        del group[name]
        group.create_dataset(name,data=getattr(matrix,name))
    group.attrs["shape"]=matrix.shape


@pytest.mark.parametrize("slip",[.25,.5,1.])
def test_independent_asymmetric_channel_signs(slip):
    c=asdict(SimulationConfig(nx=8,nz=12,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2"))
    result=channel_reference(dict(c,left_slip_length_m=slip))
    assert result["passed"] and result["relative_error"]<1e-12
    assert abs(result["left_Robin_residual_m_s"])<1e-15
    assert result["right_Dirichlet_value_m_s"]==0.


def test_real_asymmetric_native_positive_and_fixed_right(tmp_path):
    native=make_case(tmp_path/"valid_asymmetric")
    report=verify_native(native)
    assert report["passed"],report
    assert report["endpoint_trace_rows_nonzeros"]==[1,0]
    assert report["full_wall_support"]=="left_only"
    assert report["maxima"]["right_velocity_max"]==0.
    assert all(value is None or value==0. for value in comparison_metrics(native,native).values())
    with h5py.File(native,"r") as h:
        initial=h["accepted/0/eta"][-1]
        assert all(g["eta"][-1]==g["H"][1]==initial for g in h["accepted"].values())


def test_coordinate_roundoff_bound_rejects_shape_nonfinite_and_real_shift():
    expected=np.array([[-.5,.5],[-8.,-2.]])
    one_ulp=np.nextafter(expected,np.inf)
    assert _coordinate_check(one_ulp,expected,[[1.],[10.]])["passed"]
    assert not _coordinate_check(expected[:1],expected,[[1.],[10.]])["passed"]
    nonfinite=expected.copy(); nonfinite[0,0]=np.nan
    assert not _coordinate_check(nonfinite,expected,[[1.],[10.]])["passed"]
    shifted=expected.copy(); shifted[0,0]+=1e-8
    assert not _coordinate_check(shifted,expected,[[1.],[10.]])["passed"]


def test_real_native_coordinate_roundoff_is_not_a_different_mesh(tmp_path):
    native=make_case(tmp_path/"one_ulp_coordinates")
    with h5py.File(native,"r+") as h:
        h["mesh/surface_x"][3]=np.nextafter(h["mesh/surface_x"][3],np.inf)
        points=h["mesh/points"][:]
        index=int(np.flatnonzero((abs(points[0])<1.) & (points[1]<0.) & (points[1]>-10.))[0])
        h["mesh/points"][0,index]=np.nextafter(points[0,index],np.inf)
    report=verify_native(native)
    assert report["passed"],report
    assert report["mesh_coordinate_roundoff"]["surface"]["max_absolute_difference_m"]>0.
    assert report["mesh_coordinate_roundoff"]["topology_exact"]


@pytest.mark.parametrize("mutation,expected",[
    ("surface_coordinate","native_surface_mesh_mismatch"),
    ("bulk_coordinate","native_bulk_mesh_mismatch"),
    ("bulk_topology","native_bulk_mesh_mismatch"),
])
def test_real_native_mesh_change_is_not_roundoff(tmp_path,mutation,expected):
    native=make_case(tmp_path/mutation)
    with h5py.File(native,"r+") as h:
        if mutation=="surface_coordinate":
            h["mesh/surface_x"][3]+=1e-8
        elif mutation=="bulk_coordinate":
            h["mesh/points"][0,3]+=1e-8
        else:
            triangle=h["mesh/triangles"][:,0]
            h["mesh/triangles"][:,0]=triangle[[1,0,2]]
    report=verify_native(native)
    assert not report["passed"] and not report["numerical_checks_passed"],report
    assert expected in report["failures"],report


@pytest.mark.parametrize("mutation,expected",[
    ("right_vertical_unconstrained","wrong_asymmetric_free_DOFs_or_right_vertical_unconstrained"),
    ("extra_RIGHT_wall_term","forbidden_RIGHT_Navier_term_even_if_constrained"),
    ("absent_LEFT_wall_term","full_wall_operator_is_not_positive_LEFT_ONLY"),
    ("wrong_LEFT_wall_sign","full_wall_operator_is_not_positive_LEFT_ONLY"),
    ("renderer_P2_fixed_native_RIGHT_moves","RIGHT_native_contact_moves_despite_P2_label"),
    ("RIGHT_H_not_z2","RIGHT_H_not_initial_z2"),
    ("symmetric_identity","not_explicit_asymmetric_source_identity"),
    ("symmetric_contact_method","contact_is_not_asymmetric_actual_FE_trace"),
])
def test_asymmetric_native_mutations_rejected(tmp_path,mutation,expected):
    native=make_case(tmp_path/mutation)
    with h5py.File(native,"r+") as h:
        if mutation=="right_vertical_unconstrained":
            c=json.loads(h.attrs["config"]); op=_operators(c)
            right_vertical=np.intersect1d(op["fem"].wall_dofs["right"],op["fem"].component_dofs[1])
            values=h["mesh/free_velocity_dofs"][:]
            values[-1]=right_vertical[-1]
            h["mesh/free_velocity_dofs"][:]=values
        elif mutation=="extra_RIGHT_wall_term":
            c=json.loads(h.attrs["config"]); op=_operators(c); f=op["fem"]
            right=asm(scalar_mass,f.scalar.boundary("right",intorder=6)).tocoo()
            vertical=f.component_dofs[1]
            # Remove quadrature roundoff on analytically zero, non-wall traces
            # (about 1e-33 in the reduced matrix), retaining exact wall support.
            support=np.isin(vertical[right.row],f.wall_dofs["right"]) & np.isin(vertical[right.col],f.wall_dofs["right"])
            extra=coo_matrix((c["nu"]/c["left_slip_length_m"]*right.data[support],
                (vertical[right.row[support]],vertical[right.col[support]])),shape=(f.velocity.N,f.velocity.N)).tocsr()
            # This corruption changes NO reduced matrix/trajectory. A purely
            # dynamical check would miss the explicitly forbidden right law.
            reduced=extra[op["free"]][:,op["free"]]
            assert reduced.nnz==0 or np.max(abs(reduced.data))==0.
            replace_wall(h,_read_wall_matrix(h)+extra)
        elif mutation=="absent_LEFT_wall_term":
            replace_wall(h,_read_wall_matrix(h)*0.)
        elif mutation=="wrong_LEFT_wall_sign":
            replace_wall(h,-_read_wall_matrix(h))
        elif mutation=="renderer_P2_fixed_native_RIGHT_moves":
            # P2 marker metadata stays fixed, but the underlying contact moves.
            h["accepted/4/eta"][-1]+=1e-4
        elif mutation=="RIGHT_H_not_z2":
            h["accepted/4/H"][1]+=1e-4
        elif mutation=="symmetric_identity":
            ident=json.loads(h.attrs["identity"])
            ident.update(model_id="PW1_SYMMETRIC_NAVIER",source_kind="declared_Navier_linear_FEM",
                physical_model_label="DECLARED_EXTENSION — конечное скольжение Навье")
            h.attrs["identity"]=json.dumps(ident)
        elif mutation=="symmetric_contact_method":
            h.attrs["contact_method"]="actual_free_endpoint_trace_Navier"
    report=verify_native(native)
    assert not report["passed"]
    assert expected in report["failures"],report


def test_wrong_LEFT_dynamics_cannot_hide_behind_positive_saved_matrix(tmp_path,monkeypatch):
    class WrongLeftDynamics(LeftNavierFEM):
        def __init__(self,c):
            super().__init__(c)
            # Leave stored full K positive while deliberately advancing with
            # negative reduced wall friction. The stage audit must notice.
            self.K_wall=-self.K_wall
            self.K=self.K_bulk+self.K_wall
    monkeypatch.setattr(asymmetric_run,"LeftNavierFEM",WrongLeftDynamics)
    monkeypatch.setattr(LeftNavierDiagnostics,"validate",lambda *args:None)
    native=make_case(tmp_path/"wrong_left_dynamics")
    report=verify_native(native)
    assert not report["passed"]
    assert "limit:stage_momentum_relative" in report["failures"],report
    assert "limit:split_energy_relative" in report["failures"],report


def test_symmetric_native_or_mixed_left_slip_cannot_be_refinement(tmp_path):
    first=make_case(tmp_path/"left_b_half")
    second=make_case(tmp_path/"left_b_one",slip=1.)
    with pytest.raises(ValueError,match="mixed_asymmetric_physics"):
        comparison_metrics(first,second)
    assert comparison_metrics(first,second,same_physics=False)["normalized_macro_Linf"]>0.
    with h5py.File(second,"r+") as h:
        ident=json.loads(h.attrs["identity"])
        ident["model_id"]="PW1_SYMMETRIC_NAVIER"
        h.attrs["identity"]=json.dumps(ident)
    with pytest.raises(ValueError,match="symmetric_or_undeclared_native"):
        comparison_metrics(first,second,same_physics=False)


def test_real_pinned_RIGHT_slope_cannot_be_qualified_by_small_equation_residual(tmp_path):
    """Small genuine linear trajectory; a screen failure is not NOT_RUN."""
    native=make_case(tmp_path/"real_pinned_right_screen",t_end=1.)
    report=verify_native(native)
    assert report["numerical_checks_passed"] and report["passed_equations_and_BCs"],report
    assert report["maxima"]["max_slope"]>.30,report
    assert not report["linear_applicability_passed"] and not report["passed"]
    assert "limit:max_slope" in report["failures"]
    assert report["right_eta_and_H_initial_z2"] is True
    assert report["maximum_slope_witness"]["absolute_slope"]==report["maxima"]["max_slope"]
    assert report["first_slope_screen_failure"]["absolute_slope"]>.30


def test_actual_symmetric_trajectory_cannot_be_corrected_refinement(tmp_path):
    from sloshing.pinned_wetting.extension_run import execute_case as symmetric_execute, physical_identity as symmetric_physics
    from sloshing.pinned_wetting.navier import LABEL as SYMMETRIC_LABEL, NavierConfig
    directory=tmp_path/"historical_symmetric_fixture"
    directory.mkdir()
    c=NavierConfig(nx=8,nz=12,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2",slip_length_m=.5)
    identity={"run_id":"actual_symmetric_fixture","physical_model_label":SYMMETRIC_LABEL,
              "source_kind":"declared_Navier_linear_FEM","physical":symmetric_physics(c),"synthetic":False}
    symmetric_execute(ROOT,directory,c,identity,{},lambda *args:None)
    symmetric=directory/"native.h5"
    actual=make_case(tmp_path/"actual_asymmetric_fixture")
    with pytest.raises(ValueError,match="symmetric_or_undeclared_native"):
        comparison_metrics(actual,symmetric)
    report=verify_native(symmetric)
    assert not report["passed"] and not report["numerical_checks_passed"]
    assert any("forbidden_symmetric_or_right_slip_parameter" in reason for reason in report["failures"])
