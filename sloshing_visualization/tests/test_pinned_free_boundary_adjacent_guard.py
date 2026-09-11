"""Production adjacent-loop guard: actual missed loop and independent fixtures."""
import numpy as np
import pytest
from sloshing.pinned_wetting.free_boundary import adjacent_quadratic_loop, surface_has_resolved_crossing

# Five native P2 nodes of edges50/51 from the preserved failed local3 run.
# These are geometry regressions, not synthetic CFD trajectories.
GOOD=np.array([[.9922486959678646,.9921918506814649,.9923168089831718,.992676240391782,.9940109134959332],
               [.00416019893057635,.0034600336452574017,.002799301606986289,.003817492133634785,.006033874905117646]])
BAD=np.array([[.9922837820378541,.9922173621779612,.9923331727331896,.9926786655669408,.9940059585596006],
              [.004213082712039788,.003510126488850882,.0028439206598981364,.003862065247295636,.006067228559561664]])

@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('scale',[1.,1e-6,1e6])
def test_actual_accepted294_loop_cannot_hide_in_shared_chord_pair(reverse,scale):
    bad=BAD[:,::-1] if reverse else BAD
    good=GOOD[:,::-1] if reverse else GOOD
    assert adjacent_quadratic_loop(bad*scale)
    assert surface_has_resolved_crossing(bad*scale,subdivisions=2)
    assert not adjacent_quadratic_loop(good*scale)

@pytest.mark.parametrize('points',[
    [[-1.,-.5,0.,.5,1.],[0.,0.,0.,0.,0.]],
    [[-1.,-.5,0.,0.,0.],[0.,0.,0.,.5,1.]],
    [[-1.,-.5,0.,.5,1.],[1.,.25,0.,.25,1.]],
    [[-1.,-.5,0.,.5,1.],[0.,.1,.2,.1,0.]],
])
def test_only_known_join_of_straight_corner_or_curved_edges_is_allowed(points):
    assert not adjacent_quadratic_loop(np.array(points))

@pytest.mark.parametrize('points',[
    [[1.,.5,0.,.5,1.],[0.,0.,0.,0.,0.]],
    [[1.,.5,0.,.5,1.],[1.,.25,0.,.25,1.]],
    [[0.,-.5,0.,.5,1.],[0.,0.,0.,0.,0.]],
])
def test_collinear_overlap_coincident_curve_or_return_to_join_is_rejected(points):
    assert adjacent_quadratic_loop(np.array(points))
