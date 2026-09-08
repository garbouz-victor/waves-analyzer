"""Pure administrative identity, schedule-aware cost and fail-closed bindings."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pytest

SCRIPTS=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(SCRIPTS))
import step3a10_admin as admin
import step3a10_cost_audit as cost


def sample(step,core=1.,archive=.1,field=False,checkpoint=False,audit=False):
    return {"step":step,"core_s":core,"archival_s":archive,"total_s":core+archive,
            "core_class":"CORE_AUDIT" if audit else "CORE_NORMAL",
            "archive_class":cost.event_class(field,checkpoint,audit)}


def test_lineage_diff_is_separate_and_never_mutated():
    a={"execution_base_HEAD":"historical","source":"identical","stack":{"version":"same"}}
    b=copy.deepcopy(a);b["execution_base_HEAD"]="actual-current";old=copy.deepcopy(a)
    admin.require_numerical_match(a,b)
    assert a==old and b["execution_base_HEAD"]=="actual-current"


@pytest.mark.parametrize("field",["source","stack","unexpected"])
def test_any_numerical_difference_fails(field):
    a={"execution_base_HEAD":"old","source":"same","stack":"same"};b=dict(a)
    b[field]="changed"
    with pytest.raises(ValueError,match="Numerical"):admin.require_numerical_match(a,b)


def test_event_class_and_exact_future_counts():
    ev=[{"step":i,"archive_class":cost.event_class(i==2,i==5,False)} for i in range(1,7)]
    table=cost.counts(ev,3)
    assert table=={"F0_C0_A0":{"total":4,"observed":2,"future":2},
                   "F0_C1_A0":{"total":1,"observed":0,"future":1},
                   "F1_C0_A0":{"total":1,"observed":1,"future":0}}


def test_n_ge5_linear_p95():
    assert cost.estimate([1,2,3,4,10])["upper_s"]==pytest.approx(1.25*np.percentile([1,2,3,4,10],95))


@pytest.mark.parametrize("a",[[2,3],[1,4,2],[4,3,2,1]])
def test_n2to4_max_rule(a):assert cost.estimate(a)["upper_s"]==1.5*max(a)


def test_singleton_and_unseen():
    assert cost.estimate([3])["upper_s"]==6
    with pytest.raises(cost.CostAuditFailure):cost.estimate([])
    assert cost.estimate([],calibration=[1,2,3,4,5])["upper_s"]==7.5
    with pytest.raises(cost.CostAuditFailure):cost.estimate([],calibration=[1,2,3,4])


def test_unseen_core_never_uses_archive_calibration():
    train=[sample(1)];target=[sample(2,audit=True)]
    with pytest.raises(cost.CostAuditFailure):cost.predict(train,target,{"CORE_AUDIT":[1]*5,"F0_C0_A1":[1]*5})


def test_no_whole_step_p95_path_and_sparse_events_never_dropped():
    train=[sample(i,archive=.1) for i in range(1,10)]+[sample(10,archive=10,checkpoint=True)]
    target=[sample(i) for i in range(11,110)]+[sample(110,checkpoint=True)]
    p=cost.predict(train,target,{})
    expected=100*1.25+99*.125+20.
    assert p["predicted_interval_s"]==pytest.approx(expected)
    old=1.25*np.percentile([r["total_s"] for r in train],95)*len(target)
    assert old>expected
    huge=[sample(1,core=20,archive=1000,field=True,checkpoint=True,audit=True)]
    assert cost.predict(huge,[huge[0]]*3,{})["predicted_interval_s"]==3*(40+2000)


def test_backtest_underprediction_denies_authorization():
    p=cost.backtest([sample(i) for i in range(1,6)],[sample(6)],{},100.)
    assert not p["passed"]
    assert not cost.authorization_gate(100,100,[p]*4)


@pytest.mark.parametrize("full,total,expected",[(14400,19800,True),(14400.001,19000,False),
    (13000,19800.001,False),(float("nan"),100,False)])
def test_both_caps_without_slack(full,total,expected):
    assert cost.authorization_gate(full,total,[{"passed":True}]*4)==expected


def test_initialization_session_overhead_and_negative_accounting():
    obs=[sample(1),sample(2)]
    sessions=[{"status":"complete","full_path":True,"target_steps":2,"wall_s":15.,"initialization_JIT_s":10.}]
    a=cost.session_allowances(sessions,obs)
    assert a["C_init_s"]==12.5 and a["C_session_s"]==pytest.approx(4.2) and a["W_prefix_actual_s"]==15.
    sessions[0]["wall_s"]=1.
    with pytest.raises(cost.CostAuditFailure):cost.session_allowances(sessions,obs)
    with pytest.raises(cost.CostAuditFailure):cost.session_allowances([],obs)


def test_negative_core_roundoff_is_explicit_and_bounded():
    ev=[sample(1)]
    out=cost.decompose([{"step":1,"timing":{"total_s":1.,"archival_s":1.+1e-10}}],ev)
    assert out[0]["core_s"]==0 and out[0]["core_roundoff_clamped_s"]<0
    with pytest.raises(cost.CostAuditFailure):cost.decompose([{"step":1,"timing":{"total_s":1.,"archival_s":2.}}],ev)


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value))


def bound_authorization(root):
    policy={"full_L0_cap_s":14400.,"total_cap_s":19800.};write(root/"cost_audit/policy.json",policy)
    write(root/"inherited_prefix_manifest.json",{"prefix":127})
    auth={"authorized":True,"cost_policy_sha256":admin.sha(root/"cost_audit/policy.json"),
          "prefix_manifest_sha256":admin.sha(root/"inherited_prefix_manifest.json"),"schedule_sha256":"same",
          "administrative_files_sha256":{}}
    decision={"authorized":True,"policy_sha256":auth["cost_policy_sha256"],
              "prefix_manifest_sha256":auth["prefix_manifest_sha256"],"schedule_sha256":"same",
              "forecast":{"projected_complete_L0_s":100.,"projected_total_scientific_s":200.},
              "backtests":[{"passed":True}]*4}
    write(root/"cost_audit/decision.json",decision);auth["cost_decision_sha256"]=admin.sha(root/"cost_audit/decision.json")
    write(root/"continuation_authorization.json",auth)
    return auth


@pytest.mark.parametrize("file",["cost_audit/decision.json","cost_audit/policy.json","inherited_prefix_manifest.json"])
def test_hash_binding_blocks_controller(file,tmp_path):
    bound_authorization(tmp_path);admin.require_authorization(tmp_path)
    (tmp_path/file).write_text("{}")
    with pytest.raises(ValueError,match="binding"):admin.require_authorization(tmp_path)


def test_no_authorization_no_pde_import(tmp_path):
    write(tmp_path/"continuation_authorization.json",{"authorized":False})
    with pytest.raises(ValueError,match="not authorized"):admin.require_authorization(tmp_path)
    source=(SCRIPTS/"step3a10_continue.py").read_text()
    assert source.index("auth = require_authorization()")<source.index("from sloshing.multiphase.benchmarks")
    assert "d.source_guard =" not in source and "os.environ[" not in source
