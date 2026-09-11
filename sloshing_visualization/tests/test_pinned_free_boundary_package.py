"""Release guards are not numerical validation; mutations must still fail."""
import copy
import json
import numpy as np
import pytest
from sloshing.pinned_wetting.free_boundary_package import frame_map_matches_native
from sloshing.pinned_wetting.free_boundary_output import curve_bounds, tessellation_bound


@pytest.mark.parametrize("mutation",["wrong_step","wrong_time","reordered","missing","bool_step","fake5s"])
def test_frame_map_cannot_relabel_another_trajectory(mutation):
    times=np.arange(1001)*.005
    good=[{"frame":i,"time_s":i*.01,"accepted_step":2*i} for i in range(501)]
    assert frame_map_matches_native(good,times)
    bad=copy.deepcopy(good)
    if mutation=="wrong_step":bad[12]["accepted_step"]+=1
    elif mutation=="wrong_time":bad[12]["time_s"]+=.001
    elif mutation=="reordered":bad[12],bad[13]=bad[13],bad[12]
    elif mutation=="missing":bad.pop()
    elif mutation=="bool_step":bad[0]["accepted_step"]=False
    elif mutation=="fake5s":times=times*.2
    assert not frame_map_matches_native(bad,times)


def test_plot_extrema_include_curved_edge_not_only_vertices():
    # Explicit geometry fixture, not data for the mission.
    surface=np.array([[[-1.,0.,1.],[0.,1.,1.]]])
    lower,upper=curve_bounds(surface)
    assert lower[1]==0 and upper[1]==1.125
    assert tessellation_bound(surface,32)==1/2048


def test_complete_release_tests_include_validation_marker(tmp_path,monkeypatch):
    from sloshing.pinned_wetting import free_boundary_package as package
    seen=[]
    class Finished:
        returncode=0
        def poll(self):return 0
    def start(command,**kwargs):
        seen.extend(command)
        (tmp_path/"tests.xml").write_text('<testsuite tests="1" failures="0" errors="0"/>')
        return Finished()
    monkeypatch.setattr(package.subprocess,"Popen",start)
    monkeypatch.setattr(package,"tested_source_hashes",lambda root:{"fixture":"test_only"})
    monkeypatch.setattr(package,"doctor",lambda root:{"execution_HEAD":"test_only"})
    result=package.run_complete_tests(tmp_path,tmp_path,lambda:None)
    assert seen[seen.index("-o")+1]=="addopts=-ra"
    assert "../mission/pinned_wetting/solve_to_animation_v2/tests" in seen
    assert result["exit_code"]==0 and result["full_project_scope"] is True


@pytest.mark.parametrize("mutation",[None,"drop_branch_point","linearize_branch","sort_by_x","wrong_time"])
def test_export_must_preserve_entire_native_parametric_branch(tmp_path,mutation):
    import csv
    import h5py
    from sloshing.pinned_wetting.free_boundary_package import verify_full_curve_export
    from sloshing.pinned_wetting.free_boundary_output import curve
    # Synthetic overhanging geometry fixture, not solver/native mission data.
    geometry=np.array([[-1.,-.5,.8,.92,.97,.96,1.],[-.03,0.,.01,-.025,-.02,.03,.0349]])
    native=tmp_path/"test_only.h5"
    with h5py.File(native,"w") as h:
        h.create_dataset("topology/interface_nodes",data=np.arange(7))
        state=h.create_group("states/00000000")
        state.attrs["time_s"]=0.
        state.create_dataset("geometry",data=geometry)
    points=curve(geometry,32).T
    if mutation=="drop_branch_point":points=np.delete(points,-12,axis=0)
    elif mutation=="linearize_branch":points[-33:]=np.linspace(points[-33],points[-1],33)
    elif mutation=="sort_by_x":points=points[np.argsort(points[:,0])]
    exported=tmp_path/"interface_samples.csv"
    with exported.open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["time_s","vertex_id","x_m","z_m"])
        writer.writerows((.1 if mutation=="wrong_time" else 0.,i,*point) for i,point in enumerate(points))
    report=verify_full_curve_export(native,exported)
    assert report["passed"] is (mutation is None),report


@pytest.mark.parametrize("mutation",[None,"bad_exit","physics","source","looser_dt","gamma",
                                    "not_finer","no_step","stale_digest","audit_failure"])
def test_further_command_requires_compatible_finer_native_smoke(tmp_path,monkeypatch,mutation):
    import h5py
    from sloshing.pinned_wetting import free_boundary_package as package
    # Header-only synthetic fixtures; the audit is an explicit test double.
    # These files can never be used as real mission numerical evidence.
    identity={"model_id":"test_only","physical_contract":{"b":.5},
        "source_hashes":{"fixture":"synthetic"},"libraries":{"fixture":"none"},
        "method":"test_only","controls":{"nx":12,"nz":24,"dt":.005,"t_end":5.,"intorder":8,"grad_div_gamma_m2_s":.1}}
    other=copy.deepcopy(identity);other["controls"].update(nx=24,nz=48)
    if mutation=="physics":other["physical_contract"]["b"]=1.
    if mutation=="source":other["source_hashes"]["fixture"]="wrong"
    if mutation=="looser_dt":other["controls"]["dt"]=.01
    if mutation=="gamma":other["controls"]["grad_div_gamma_m2_s"]=.2
    original,refined=tmp_path/"original_test_only.h5",tmp_path/"refined_test_only.h5"
    for path,ident,count in ((original,identity,2),(refined,other,2 if mutation=="not_finer" else 4)):
        with h5py.File(path,"w") as h:
            h.attrs["numerical_identity"]=json.dumps(ident)
            h.attrs["last_accepted_step"]=0 if mutation=="no_step" else 1
            h.create_dataset("topology/triangles",data=np.zeros((3,count),dtype=int))
    record={"native":str(refined),"native_sha256":package.sha256(refined),
        "command":["test_only","run","--free-boundary","--stop-after","1"],
        "exit_code":1 if mutation=="bad_exit" else 0}
    if mutation=="stale_digest":record["native_sha256"]="wrong"
    package.write_json(tmp_path/"further_refinement_command.json",record)
    monkeypatch.setattr(package,"verify_native",lambda *args,**kwargs:{"passed":mutation!="audit_failure"})
    answer=package.further_refinement_evidence(tmp_path,tmp_path,original,lambda:None)
    assert answer["passed"] is (mutation is None),answer


def test_no_unrun_further_refinement_command_can_pass(tmp_path):
    from sloshing.pinned_wetting.free_boundary_package import further_refinement_evidence
    assert not further_refinement_evidence(tmp_path,tmp_path,tmp_path/"unused",lambda:None)["passed"]
