import json
from pathlib import Path
import numpy as np
import pytest
from sloshing.multiphase.linearized_ch import json_hash
from sloshing.multiphase.phase_rate_schedule import RateSchedule, cost_projection, scalar_convergence
from sloshing.multiphase.phase_rate_history import (Journal,archive,load_archive,atomic_json,recover_prefix,
    WriterLock,validate_complete)


def tiny_plan(dt=1e-6, counts=(3,3)):
    blocks=[]; start=0.
    for i,n in enumerate(counts):
        h=dt*2**i; end=start+n*h
        blocks.append(dict(t_start=start,t_end=end,dt=h,steps=n)); start=end
    p=dict(blocks=blocks,t_end=start,total_steps=sum(counts))
    return {**p,"plan_sha256":json_hash(p)}


def test_historical_schedule_counts_hashes_and_exact_nested_endpoints():
    plan=json.loads(Path("validation_results/step3a4/schedule_optimization/optimized_plan.json").read_text())
    sha="5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8"
    levels=[RateSchedule(plan,f,expected_hash=sha,expected_count=1068) for f in (1,2,4)]
    assert [s.nsteps for s in levels]==[1068,2136,4272]
    assert levels[0].block_ends[0]==122 and levels[0].pilot_end==127
    ends=[r.end for r in levels[0].intervals]
    for s in levels:
        assert s.sha256==RateSchedule(plan,s.factor).sha256
        assert [r.end for r in s.intervals[s.factor-1::s.factor]]==ends
        assert all(a.end==b.start for a,b in zip(s.intervals,s.intervals[1:]))
        assert s.intervals[-1].end==1e-4
        assert [s.intervals[i*s.factor].dt for i in range(1068)]==[r.dt/s.factor for r in levels[0].intervals]
        assert s.common_parent_indices==levels[0].common_parent_indices
    with pytest.raises(ValueError): RateSchedule(plan,expected_hash="wrong")
    with pytest.raises(ValueError): RateSchedule(plan,expected_count=898)


def test_cost_never_uses_old_full_CHNS_and_does_not_relax_budget():
    limits={"isolated_series_wall_s":100.}
    rows=[dict(step=i,timing={"total_s":1.}) for i in range(1,5)]
    p=cost_projection(rows,[10,20,30],5.,1.,limits)
    assert p["projected_s_per_step"]==1.25 and p["series_cost_ok"]
    assert not cost_projection(rows,[100],5.,1.,limits)["series_cost_ok"]
    with pytest.raises(ValueError): cost_projection([], [1],0,0,limits)


def test_scalar_BE_refinement_and_observable_failure():
    I=[]; dF=[]; B=[]; fields=[]
    for f in (1,2,4):
        q=1.; integ=0.
        for r in RateSchedule(tiny_plan(.05,(3,3)),f).intervals:
            new=q/(1+r.dt); integ+=r.dt/2*(q*q+new*new); q=new
        fields.append(q); I.append(integ); dF.append(.5*(q*q-1)); B.append(dF[-1]+integ)
    gates={"last_D_difference_relative":.05,"final_budget_relative":.05}
    gaps=abs(np.diff(fields)); energies=np.array(dF)+.5
    result=scalar_convergence(I,dF,B,gaps,energies,gates)
    assert result["passed"]
    assert not scalar_convergence(I,dF,B,[1.,2.],energies,gates)["passed"]
    with pytest.raises(ValueError): scalar_convergence(I[:2],dF[:2],B[:2],gaps,energies,gates)


def test_journal_and_checkpoint_corruption_rejected(tmp_path):
    j=Journal(tmp_path/"scalar_history.jsonl")
    j.append(dict(step=0,time=0.)); j.append(dict(step=1,time=1.))
    p=j.prefix(); meta=dict(identity={"source":"a"},accepted_steps=1,journal_prefix=p,accepted=True)
    archive(tmp_path/"checkpoint",{"phase_rate":np.array([.123456789])},meta)
    arrays,loaded=load_archive(tmp_path/"checkpoint",{"source":"a"})
    assert arrays["phase_rate"][0]==.123456789
    with pytest.raises(ValueError): load_archive(tmp_path/"checkpoint",{"source":"b"})
    j.append(dict(step=2,time=2.))
    restored=recover_prefix(tmp_path,loaded)
    assert len(restored.rows)==2 and list(tmp_path.glob("orphan-*"))
    atomic_json(tmp_path/"FAILED.json",{"reason":"negative SNES"})
    with pytest.raises(ValueError,match="failure"): recover_prefix(tmp_path,loaded)
    with pytest.raises(FileNotFoundError): validate_complete(tmp_path,{},RateSchedule(tiny_plan()))


def test_atomic_archives_are_independent_and_writer_exclusive(tmp_path):
    rate=np.array([1.,2.]); archive(tmp_path/"a",{"phase_rate":rate},{"accepted":True})
    rate[:]=0.
    assert np.array_equal(load_archive(tmp_path/"a")[0]["phase_rate"],[1.,2.])
    with pytest.raises(FileExistsError): archive(tmp_path/"a",{}, {})
    a=WriterLock(tmp_path)
    try:
        with pytest.raises(RuntimeError,match="Concurrent"): WriterLock(tmp_path)
    finally: a.close()


def test_journal_does_not_alias_nested_timing_and_torn_suffix_is_preserved(tmp_path):
    path=tmp_path/"scalar_history.jsonl"; journal=Journal(path)
    row={"step":0,"timing":{"seconds":1.}}
    journal.append(row); row["timing"]["seconds"]=9.
    assert journal.rows==Journal.read(path)
    prefix=journal.prefix()
    with path.open("ab") as stream: stream.write(b'{"torn":')
    with pytest.raises(ValueError,match="Torn"): Journal.read(path)
    restored=recover_prefix(tmp_path,{"journal_prefix":prefix,"accepted_steps":0})
    assert restored.rows[0]["timing"]["seconds"]==1.
    assert list(tmp_path.glob("orphan-*/scalar_history.jsonl"))


def test_physical_fail_cannot_be_overridden_by_solver_or_comparison_pass():
    from sloshing.multiphase.phase_rate_trajectory import physical_gates,finite_tree
    row=dict(mass_error_relative_to_domain=2e-10,phase_mass=0.,wall_crossings=2,delta_E_interval=0.)
    p=dict(nonlinear_mass_domain=1e-10,certified_cells=8,energy_growth_J_per_m=1e-9,weak_work_J_per_m=1e-12)
    assert not all(physical_gates(row,{"cells_across_transition_certified_min":9},None,p,0.,1.,2).values())
    assert not finite_tree({"bad":[float("nan")]})
