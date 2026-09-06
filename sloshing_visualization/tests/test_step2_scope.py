from sloshing.visualization.scope import load_scope, region_codes, event_frames
import numpy as np


def test_scope_is_separate_and_boundaries_are_explicit():
    scope=load_scope()
    assert scope["interpretation"]=="CONDITIONALLY QUALIFIED FOR BULK MODEL VISUALIZATION"
    np.testing.assert_array_equal(region_codes([-.99,-.98,-.9,0,.9,.98,.99],scope),[2,1,0,0,0,1,2])


def test_key_frames_use_measured_events_and_saved_times():
    modal={"first_zero_crossing_s":.413,"first_opposite_extremum":{"time":.811},"next_zero_crossing_s":1.227,
           "extrema":[{"time":1.619,"kind":"maximum"},{"time":4.043,"kind":"minimum"}]}
    rows=event_frames(np.arange(1001)*.005,modal)
    assert rows[1]["selected_time_s"]==.41500000000000004 or np.isclose(rows[1]["selected_time_s"],.415)
    assert rows[-1]["event_time_s"]==4.043


def test_storyboard_includes_several_measured_late_extrema():
    from sloshing.visualization.storyboard import storyboard_events
    times=np.arange(1001)*.005
    meta={"key_frames":[{"index":0,"name":"t0"}],"modal":{"extrema":[{"time":t,"kind":"maximum"} for t in (.8,2.4,3.2,4.8)]}}
    rows=storyboard_events(times,meta)
    assert len(rows)==4 and [r["index"] for r in rows]==[0,480,640,960]
