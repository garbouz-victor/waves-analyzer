"""PW2 implementation checks; synthetic manufactured states are tests only."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import h5py
import pytest
from sloshing.pinned_wetting import free_boundary as fb
from sloshing.pinned_wetting import free_boundary_run as runner

ROOT=Path(__file__).resolve().parents[2]


def test_initial_line_full_trace_and_asymmetric_operators():
    c=fb.Controls(nx=4,nz=8);m=fb.initial_mesh(c);o=fb.Operators(m,c)
    ids=fb.interface_nodes(m)
    np.testing.assert_allclose(m.p[1,ids],m.p[0,ids]*np.tan(np.deg2rad(2)),atol=1e-14)
    assert o.components[1][ids[0]] in o.free
    assert o.components[1][ids[-1]] in o.fixed
    assert np.max(abs(o.Kl[o.wall["right"]].data),initial=0.)==0.
    assert not np.intersect1d(o.free,o.wall["right"]).size
    assert not np.intersect1d(o.free,o.wall["bottom"]).size
    assert o.Kl.diagonal().min()>=0.


def test_grad_div_full_form_psd_and_separate_signed_work():
    from skfem import BilinearForm,asm
    from skfem.helpers import div
    c=fb.Controls(nx=4,nz=8,hydrostatic_split=True,skew_divergence=True,grad_div_gamma_m2_s=.1)
    m=fb.initial_mesh(c);o=fb.Operators(m,c)
    @BilinearForm
    def independent_grad_div(u,v,w): return .1*div(u)*div(v)
    reference=asm(independent_grad_div,o.vbasis)
    np.testing.assert_allclose(o.Kd.toarray(),reference.toarray(),rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(o.Kd.toarray(),o.Kd.T.toarray(),rtol=0,atol=1e-13)
    rng=np.random.default_rng(604)
    for _ in range(5):
        v=rng.normal(size=o.vbasis.N)
        assert v@(o.Kd@v)>=0
    zero=fb.Operators(m,replace(c,grad_div_gamma_m2_s=0.))
    assert zero.Kd.nnz==0
    np.testing.assert_array_equal(o.fixed,zero.fixed)
    np.testing.assert_array_equal(o.Kl.toarray(),zero.Kl.toarray())
    _,orientation=fb.quadratic_minima(m)
    m1,v1,p,row=fb.step(m,np.zeros(o.vbasis.N),c.dt,c,orientation)
    mid=fb.Operators(replace(m,doflocs=(m.p+m1.p)*.5),c,v1*.5)
    direct=-c.dt*float((v1*.5)@(mid.Kd@(v1*.5)))
    assert direct<0
    assert row["numerical_grad_div_work"]==pytest.approx(direct,rel=1e-13)
    assert row["bulk_loss"]>=0 and row["left_wall_loss"]>=0
    assert np.max(abs(v1[o.fixed]))==0


@pytest.mark.parametrize("gamma",[-.1,float("nan"),float("inf")])
def test_grad_div_coefficient_rejects_invalid_numerical_values(gamma):
    with pytest.raises(ValueError): fb.Controls(grad_div_gamma_m2_s=gamma)


def test_flat_hydrostatic_rest_is_solved_not_prescribed():
    c=fb.Controls(nx=4,nz=8);m=fb.initial_mesh(c,0.);o=fb.Operators(m,c)
    a,p=o.startup();assert max(abs(a))<1e-10
    _,sign=fb.quadratic_minima(m)
    m1,v,p,row=fb.step(m,np.zeros(o.vbasis.N),c.dt,c,sign)
    assert max(abs(v))<1e-11
    np.testing.assert_allclose(m1.p,m.p,rtol=0,atol=2e-14)
    assert abs(row["volume_m2"]-20.)<1e-12


def test_full_material_kinematics_and_mass():
    c=fb.Controls(nx=4,nz=8);m=fb.initial_mesh(c);o=fb.Operators(m,c)
    _,sign=fb.quadratic_minima(m);v0=np.zeros(o.vbasis.N)
    m1,v1,p,row=fb.step(m,v0,c.dt,c,sign)
    np.testing.assert_allclose(m1.p-m.p,c.dt*o.nodal_vector(v1+v0)/2,atol=1e-15)
    assert abs(row["volume_m2"]-20.)<1e-12
    assert row["bulk_loss"]>=0 and row["left_wall_loss"]>=0
    assert row["momentum_residual"]<1e-9
    assert np.max(abs(v1[o.wall["right"]]))==0
    ids=fb.interface_nodes(m)
    assert m1.p[1,ids[-1]]==m.p[1,ids[-1]]
    assert m1.p[1,ids[0]]>m.p[1,ids[0]]


def test_exact_jacobian_interior_and_gcl_translation():
    c=fb.Controls(nx=4,nz=8);m=fb.initial_mesh(c)
    minimum,orientation=fb.quadratic_minima(m)
    shifted=replace(m,doflocs=m.p+np.array([[.21],[-.13]]))
    shifted_min,shifted_sign=fb.quadratic_minima(shifted)
    np.testing.assert_allclose(minimum,shifted_min,rtol=1e-12,atol=1e-14)
    np.testing.assert_array_equal(orientation,shifted_sign)
    # A P2 mid-edge fold with all three vertices intact is not a valid cell.
    bad=m.p.copy();s=fb.Operators(m,c).scalar
    element=s.element_dofs[:,0];bad[:,element[3]]=bad[:,element[0]]+5*(bad[:,element[2]]-bad[:,element[0]])
    folded,_=fb.quadratic_minima(replace(m,doflocs=bad))
    assert folded.min()<0


def test_checkpoint_and_HEAD_are_separate(tmp_path,monkeypatch):
    monkeypatch.setattr(runner,"REL_RUNS",str(tmp_path/"runs"))
    actual=runner.doctor(ROOT);head=["creation_HEAD"]
    monkeypatch.setattr(runner,"doctor",lambda root:{**actual,"execution_HEAD":head[0]})
    c=fb.Controls(nx=4,nz=8,dt=.0025,t_end=.0075,grad_div_gamma_m2_s=.1)
    state={};noop=lambda *a:None
    partial=runner.run_case(ROOT,c,"test",state,noop,stop_after=1)
    assert partial["accepted_step"]==1
    with h5py.File(partial["native"]) as h:before=h["states/00000001/geometry"][:]
    head[0]="documentation_only_new_HEAD"
    complete=runner.run_case(ROOT,c,"test",state,noop)
    assert complete["completed"] and len(state["run_ids"])==1
    with h5py.File(complete["native"]) as h:
        np.testing.assert_array_equal(h["states/00000001/geometry"][:],before)
        coating=json.loads(h["states/00000003"].attrs["coating"])
        history=json.loads(h["states/00000003"].attrs["diagnostics"])
        assert history["cumulative_numerical_grad_div_work"]<0
        assert coating["accepted_step"]==3 and coating["H"][0]>=partial["coating"]["H"][0]
        assert coating["markers"][0]==partial["coating"]["markers"][0]
    provenance=json.loads((Path(complete["native"]).parent/"provenance.json").read_text())
    assert [p["HEAD"] for p in provenance["execution_HEAD_history"]]==["creation_HEAD","documentation_only_new_HEAD"]
    digest=runner.sha256(complete["native"])
    reused=runner.run_case(ROOT,c,"test",state,noop)
    assert reused==complete and runner.sha256(complete["native"])==digest


def test_conservative_rezone_keeps_entire_interface_and_momentum():
    from sloshing.pinned_wetting.free_boundary_remesh import rezone
    c=fb.Controls(nx=4,nz=8,hydrostatic_split=True,skew_divergence=True)
    reference=fb.initial_mesh(c);o=fb.Operators(reference,c)
    _,sign=fb.quadratic_minima(reference)
    m,v,p,row=fb.step(reference,np.zeros(o.vbasis.N),c.dt,c,sign)
    rezoned,transferred,evidence=rezone(reference,m,v,c,sign)
    ids=fb.interface_nodes(reference)
    np.testing.assert_array_equal(m.p[:,ids],rezoned.p[:,ids])
    assert abs(evidence["volume_change"])<1e-12
    assert evidence["momentum_error"]<1e-12
    assert evidence["inverse_map_error"]<1e-10
    assert evidence["projection_residual"]<1e-9
    r=fb.Operators(rezoned,c)
    assert np.max(abs(transferred[r.fixed]))==0.
    assert np.linalg.norm(r.B@transferred,np.inf)<1e-10


def test_checkpoint_on_output_time_serializes_native_boolean(tmp_path,monkeypatch):
    monkeypatch.setattr(runner,"REL_RUNS",str(tmp_path/"runs"))
    c=fb.Controls(nx=4,nz=8,dt=.005,t_end=.02)
    result=runner.run_case(ROOT,c,"outputtime",{},lambda *a:None,stop_after=2)
    assert result["time_s"]==.01 and result["completed"] is False
    assert json.loads(json.dumps(result))["target_5s"] is False


def test_material_inertia_contains_nonzero_eulerian_convection():
    """Manufactured transport identity only, NOT a physical vessel trajectory.

    Steady Eulerian extensional v=A*x has ∂t v=0 but (v·grad)v=A²*x.
    The material coefficient derivative must reproduce this nonzero term.
    """
    from skfem import LinearForm,asm
    c=fb.Controls(nx=4,nz=8);m0=fb.initial_mesh(c)
    dt=.02;k=.4;rate=np.array([[k],[-k]])
    X1=(1+dt*rate/2)/(1-dt*rate/2)*m0.p
    mm=replace(m0,doflocs=(m0.p+X1)/2);o=fb.Operators(mm,c)
    v0=np.zeros(o.vbasis.N);v1=v0.copy()
    for j,indices in enumerate(o.components):
        v0[indices]=rate[j,0]*m0.p[j];v1[indices]=rate[j,0]*X1[j]
    @LinearForm
    def convection(v,w): return k*k*(w.x[0]*v[0]+w.x[1]*v[1])
    measured=o.M@((v1-v0)/dt)
    expected=asm(convection,o.vbasis)
    assert np.linalg.norm(expected)>1e-3
    np.testing.assert_allclose(measured,expected,rtol=1e-11,atol=2e-13)


def test_parametric_free_traction_normal_and_left_only_friction_sign():
    c=fb.Controls(nx=4,nz=8,hydrostatic_split=True)
    m=fb.initial_mesh(c);o=fb.Operators(m,c);slope=np.tan(np.deg2rad(2.))
    # Independent integrals on z=k*x, n ds=(-k,1)dx, x in[-1,1].
    for component,expected in ((0,9.81*slope*slope*2/3),(1,-9.81*slope*2/3)):
        test=np.zeros(o.vbasis.N)
        test[o.components[component]]=o.scalar.doflocs[0]
        np.testing.assert_allclose(o.f@test,expected,rtol=1e-12,atol=1e-14)
    w=np.zeros(o.vbasis.N);w[o.components[1]]=1.
    np.testing.assert_allclose(w@(o.Kl@w),.02*(10.-slope),rtol=1e-12)
    assert np.max(abs(o.Kl[o.wall["right"]].data),initial=0.)==0.


def test_curved_inverse_fallback_has_no_outside_extrapolation():
    from types import SimpleNamespace
    from sloshing.pinned_wetting.free_boundary_remesh import shape,_bounded_inverse_fallback
    # Manufactured curved triangle only; not a vessel surface or solver output.
    geometry=np.array([[0.,1.,0.,.55,.6,-.1],[0.,0.,1.,-.1,.55,.6]])
    reference=np.array([[.001,.91,.2],[.6,.03,.2]])
    points=np.column_stack((geometry@shape(reference)[0],[3.,3.]))
    cells=np.full(4,-1,dtype=int);refs=np.zeros((2,4));errors=np.zeros(4)
    _bounded_inverse_fallback(SimpleNamespace(p=geometry),np.arange(6)[:,None],points,cells,refs,errors)
    np.testing.assert_array_equal(cells,[0,0,0,-1])
    np.testing.assert_allclose(refs[:,:3],reference,rtol=0,atol=1e-9)
    np.testing.assert_allclose(geometry@shape(refs[:,:3])[0],points[:,:3],rtol=0,atol=1e-10)


def test_record_detector_cannot_beat_an_unmarked_higher_plateau():
    # Synthetic contact histories, not numerical trajectories.
    assert not runner.resolved_left_record(.02,.025,.022,.03,-.035)
    assert not runner.resolved_left_record(.02,.025,.022,np.nextafter(.025,np.inf),-.035)
    assert runner.resolved_left_record(.02,.031,.022,.031,-.035)


def test_actual_095_curved_cell_locator_regression(monkeypatch):
    from types import SimpleNamespace
    from sloshing.pinned_wetting import free_boundary_remesh as remap
    # Geometry/target extracted read-only from material-spacefine-5451846be9bc
    # at .95s, old cell1565. The affine test velocity is manufactured ONLY
    # to check interpolation; it is not that vessel's physical velocity.
    geometry=np.array([
        [.9829764897173254,.9866258613452504,1.,.9845541807150505,.9935568239823769,.9926341052982],
        [-.0255618272877031,-.02024394118464344,-.0010939153091737356,-.022786823752503727,-.01284008180411124,-.017046694508788544]])
    point=np.array([[.9966377343364423],[-.00904403977738791]])
    mesh=SimpleNamespace(p=geometry,t=np.arange(3)[:,None])
    old=SimpleNamespace(mesh=mesh,scalar=SimpleNamespace(element_dofs=np.arange(6)[:,None]),nodal_vector=lambda v:v)
    velocity=np.array([geometry[0]+2*geometry[1],3*geometry[0]-geometry[1]])
    original=remap._bounded_inverse_fallback
    monkeypatch.setattr(remap,"_bounded_inverse_fallback",lambda *args:None)
    with pytest.raises(ValueError,match="unlocated quadrature"):
        remap.locate_and_evaluate(old,velocity,point)
    monkeypatch.setattr(remap,"_bounded_inverse_fallback",original)
    values,cells,refs,error=remap.locate_and_evaluate(old,velocity,point)
    np.testing.assert_array_equal(cells,[0])
    np.testing.assert_allclose(values,np.array([point[0]+2*point[1],3*point[0]-point[1]]),atol=3e-10,rtol=0)
    assert error<1e-10 and refs.min()>=0 and refs.sum()<=1.
def test_early_crossing_guard_preserves_parameter_order_and_shared_nodes():
    from sloshing.pinned_wetting.free_boundary import surface_has_resolved_crossing
    # Synthetic piecewise straight P2 curve; nonadjacent edges cross.
    crossing=np.array([[-1.,-.4,.2,0.,-.2,0.,.2,.6,1.],
                       [0.,.1,.2,.2,.2,0.,-.2,-.05,.1]])
    assert surface_has_resolved_crossing(crossing)
    ordinary=np.array([[-1.,-.5,0.,.5,1.],[0.,.1,.2,.1,0.]])
    assert not surface_has_resolved_crossing(ordinary)
    assert not surface_has_resolved_crossing(ordinary[:,::-1])
