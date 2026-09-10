"""Corrective production constraints, coating continuation and commit provenance."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
import yaml

from sloshing.pinned_wetting.asymmetric import LeftNavierConfig, LeftNavierFEM, channel_benchmark, MODEL_ID
from sloshing.pinned_wetting.navier import NavierIntegrator
from sloshing.pinned_wetting.asymmetric_run import execute_case, run_case

ROOT = Path(__file__).resolve().parents[2]


def controls(**kw):
    return LeftNavierConfig(**dict(nx=8,nz=12,dt=.005,t_end=.02,snapshot_dt=.01,integrator="sdirk2",**kw))


@pytest.mark.parametrize("b",[.25,.5,1.])
def test_asymmetric_analytic_channel(b):
    result = channel_benchmark(controls(left_slip_length_m=b))
    assert result["relative_error"] < 1e-10
    assert result["right_velocity_max"] == 0.
    assert result["analytic_left_robin_residual"] < 1e-14


def test_left_free_right_fixed_exact_constraints():
    c=controls(); f=LeftNavierFEM(c)
    assert [f.R.getrow(i).nnz for i in (0,f.R.shape[0]-1)] == [1,0]
    assert not np.intersect1d(f.free,f.wall_dofs["right"]).size
    right=f.wall_dofs["right"]
    assert np.max(abs(f.K_wall_full[right].data),initial=0) == 0.
    solver=NavierIntegrator(f); initial=solver.initial_state(); state=initial
    for _ in range(c.nsteps):
        state=solver.advance(state)
        v=f.expand_velocity(state.velocity)
        assert state.eta[-1] == initial.eta[-1]
        assert np.max(abs(v[right])) == 0.
        assert np.max(abs(v[f.wall_dofs["bottom"]])) == 0.
    assert abs(state.eta[0]-initial.eta[0])>1e-9
    assert state.wall_dissipated_energy>0 and state.bulk_dissipated_energy>0


def test_asymmetric_zero_equilibrium():
    c=controls(alpha_deg=0.); solver=NavierIntegrator(LeftNavierFEM(c)); state=solver.initial_state()
    for _ in range(c.nsteps): state=solver.advance(state)
    assert np.max(abs(state.eta)) == np.max(abs(state.velocity)) == 0.


def test_asymmetric_continuation_native_and_coating(tmp_path):
    c=controls(); identity={"run_id":"ASYMMETRIC_CONTINUATION","model_id":MODEL_ID}
    a,b=tmp_path/"whole",tmp_path/"split"; a.mkdir(); b.mkdir(); noop=lambda *a:None
    execute_case(ROOT,a,c,identity,{},noop)
    execute_case(ROOT,b,c,identity,{},noop,stop_after=2)
    execute_case(ROOT,b,c,identity,{},noop)
    with h5py.File(a/"native.h5") as aa,h5py.File(b/"native.h5") as bb:
        for key in aa["accepted"]:
            ga,gb=aa["accepted"][key],bb["accepted"][key]
            for field in ga: np.testing.assert_array_equal(ga[field][:],gb[field][:])
            assert ga.attrs["coating"]==gb.attrs["coating"]
            assert ga["H"][1] == aa["accepted/0/eta"][-1]


def test_changed_HEAD_reuses_compatible_checkpoint_and_records_provenance(tmp_path,monkeypatch):
    import sloshing.pinned_wetting.asymmetric_run as runner
    cfg=yaml.safe_load((ROOT/"mission/pinned_wetting/left_navier_right_noslip/CONFIG.yaml").read_text())
    cfg["runs_root"]=str(tmp_path/"runs")
    head=["creation_commit"]
    monkeypatch.setattr(runner,"doctor",lambda root:{"execution_HEAD":head[0]})
    state={"run_ids":[]}; noop=lambda *a:None; c=controls()
    first=run_case(ROOT,cfg,c,"head_regression",state,noop,stop_after=2)
    directory=Path(first["native"]).parent
    identity=(directory/"identity.json").read_bytes()
    head[0]="documentation_only_commit"
    complete=run_case(ROOT,cfg,c,"head_regression",state,noop)
    assert complete["completed"] and len(state["run_ids"])==1
    assert (directory/"identity.json").read_bytes()==identity
    provenance=json.loads((directory/"provenance.json").read_text())
    assert [entry["HEAD"] for entry in provenance["execution_HEAD_history"]]==["creation_commit","documentation_only_commit"]
    digest=(directory/"native.h5").read_bytes()
    reused=run_case(ROOT,cfg,c,"head_regression",state,noop)
    assert reused==complete and (directory/"native.h5").read_bytes()==digest


def test_pyyaml_is_a_declared_runtime_dependency():
    requirements=(ROOT/"sloshing_visualization/requirements.txt").read_text()
    project=(ROOT/"sloshing_visualization/pyproject.toml").read_text()
    assert "PyYAML>=6,<7" in requirements
    assert '"PyYAML>=6,<7"' in project.split("[project.optional-dependencies]")[0]


def test_actual_cli_further_refinement_survives_real_documentation_commit(tmp_path):
    """Real CLI, real temporary git commits; never mutates the user's repository."""
    root=tmp_path/"repo"; root.mkdir()
    mission=root/"mission/pinned_wetting"
    (mission/"left_navier_right_noslip").mkdir(parents=True)
    shutil.copyfile(ROOT/"mission/pinned_wetting/CONFIG.yaml",mission/"CONFIG.yaml")
    cfg=yaml.safe_load((ROOT/"mission/pinned_wetting/left_navier_right_noslip/CONFIG.yaml").read_text())
    cfg.update(baseline_mesh=[8,12],pilot_t_end_s=.02)
    (mission/"left_navier_right_noslip/CONFIG.yaml").write_text(yaml.safe_dump(cfg))
    (mission/"state.json").write_text(json.dumps({"run_ids":[],"active_jobs":[],
        "budget":{"actual_local_job_wall_s":0.,"auxiliary_wall_time_upper_bound_s":0.}}))
    scripts=root/"sloshing_visualization/scripts"; scripts.mkdir(parents=True)
    shutil.copyfile(ROOT/"sloshing_visualization/scripts/pinned_wetting_mission.py",scripts/"pinned_wetting_mission.py")
    (root/"sloshing_visualization/src").symlink_to(ROOT/"sloshing_visualization/src",target_is_directory=True)
    def git(*args):
        return subprocess.check_output(["git","-c","user.name=PW1 test","-c","user.email=pw1@example.invalid",*args],
            cwd=root,text=True,stderr=subprocess.STDOUT).strip()
    git("init","-q"); (root/"note.txt").write_text("creation\n")
    git("add","note.txt"); git("commit","-qm","fixture creation")
    creation=git("rev-parse","HEAD")
    env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1","OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"}
    def cli(command,*extra):
        result=subprocess.run([sys.executable,str(scripts/"pinned_wetting_mission.py"),command,
            "--corrected","--mesh-level","1","--dt-scale","0.5",*extra],cwd=root,env=env,
            capture_output=True,text=True,timeout=90)
        assert result.returncode==0, result.stdout+result.stderr
    cli("pilot","--stop-after","2")
    run=next((root/cfg["runs_root"]).iterdir())
    initial_identity=(run/"identity.json").read_bytes()
    with h5py.File(run/"native.h5") as h:
        assert sorted(h["accepted"],key=int)==["0","1","2"]
        original_eta=h["accepted/2/eta"][:]
        original_H=h["accepted/2/H"][:]
    saved=json.loads((mission/"state.json").read_text())
    assert saved["accepted_step"]==2 and saved["corrected_runs"]["pilot"]["completed"] is False
    (root/"note.txt").write_text("documentation changed; equations unchanged\n")
    git("add","note.txt"); git("commit","-qm","documentation only")
    updated=git("rev-parse","HEAD"); assert updated!=creation
    cli("resume","--stop-after","4")
    with h5py.File(run/"native.h5") as h:
        assert sorted(h["accepted"],key=int)==["0","1","2","3","4"]
        np.testing.assert_array_equal(h["accepted/2/eta"][:],original_eta)
        np.testing.assert_array_equal(h["accepted/2/H"][:],original_H)
    saved=json.loads((mission/"state.json").read_text())
    assert saved["accepted_step"]==4 and saved["corrected_runs"]["pilot"]["completed"] is False
    cli("pilot")
    assert len(list((root/cfg["runs_root"]).iterdir()))==1
    assert (run/"identity.json").read_bytes()==initial_identity
    provenance=json.loads((run/"provenance.json").read_text())
    assert [v["HEAD"] for v in provenance["execution_HEAD_history"]]==[creation,updated]
    with h5py.File(run/"native.h5") as h:
        c=json.loads(h.attrs["config"])
        assert (c["nx"],c["nz"],c["dt"],c["left_slip_length_m"])==(16,24,.0025,.5)
        assert h["accepted/8"].attrs["time_s"]==.02
    before=(run/"native.h5").read_bytes(); cli("pilot")
    assert before==(run/"native.h5").read_bytes()


def test_cli_reuses_integer_float_equivalent_controls_without_duplicate(tmp_path):
    from sloshing.pinned_wetting.asymmetric_mission import pilot_case, controls as cli_controls
    cfg=yaml.safe_load((ROOT/"mission/pinned_wetting/left_navier_right_noslip/CONFIG.yaml").read_text())
    original=yaml.safe_load((ROOT/"mission/pinned_wetting/CONFIG.yaml").read_text())
    cfg.update(runs_root=str(tmp_path/"runs"),baseline_mesh=[8,12],pilot_t_end_s=1)
    state={"run_ids":[]}; noop=lambda *a:None
    c=cli_controls(original,cfg)
    first=run_case(ROOT,cfg,c,"pilot_space",state,noop,stop_after=2)
    cfg["pilot_t_end_s"]=1.0
    again=pilot_case(ROOT,original,cfg,SimpleNamespace(mesh_level=0,dt_scale=1.,stop_after=4),state,noop)
    assert first["native"]==again["native"]
    assert len(list((tmp_path/"runs").iterdir()))==1
    with h5py.File(again["native"]) as h:
        assert sorted(h["accepted"],key=int)==["0","1","2","3","4"]
