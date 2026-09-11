"""Exact candidate-pruned comparisons against the old brute distance formula.

All arrays are explicitly synthetic geometry tests, not numerical trajectories.
"""
import numpy as np
import pytest

from sloshing.pinned_wetting import free_boundary_compare as comparison


def brute_reference(points, target):
    start = target[:-1]
    edge = target[1:]-start
    denominator = np.maximum(np.sum(edge*edge, axis=1), 1e-300)
    maximum = 0.
    for first in range(0, len(points), 128):
        point = points[first:first+128, None, :]
        fraction = np.clip(np.sum((point-start)*edge, axis=2)/denominator, 0., 1.)
        difference = point-start-fraction[:, :, None]*edge
        value = float(np.sqrt(np.max(np.min(np.sum(difference**2, axis=2), axis=1))))
        maximum = max(maximum, value)
    return maximum


@pytest.mark.parametrize("seed", [11, 613, 48891])
def test_pruned_directed_distance_matches_brute_for_every_random_point(seed):
    rng = np.random.default_rng(seed)
    target = rng.normal(size=(137, 2)).cumsum(axis=0)
    points = rng.normal(size=(263, 2))*4+target.mean(axis=0)
    assert comparison._directed(points, target) == brute_reference(points, target)
    # A correct maximum cannot conceal a wrongly pruned individual query.
    for point in points[::13]:
        assert comparison._directed(point[None, :], target) == brute_reference(point[None, :], target)


@pytest.mark.parametrize("target", [np.array([[.2, -.3], [.2, -.3]]),
                                     np.repeat([[.2, -.3]], 80, axis=0),
                                     np.array([[0., 0.], [1., 0.], [1., 0.], [2., 0.]])])
def test_pruned_distance_handles_zero_length_and_repeated_segments(target):
    points = np.array([[.2, -.3], [0., 0.], [1.1, 2.], [-.1, -.4]])
    assert comparison._directed(points, target) == brute_reference(points, target)


def test_long_segment_cannot_be_pruned_by_nearby_short_segment_midpoints():
    target = np.array([[0., 0.], [1000., 0.], [1000., 10.], [.55, .3], [.6, .3]])
    points = np.array([[.5, .01], [.51, .02], [.52, .005]])
    assert comparison._directed(points, target) == brute_reference(points, target)
    assert comparison._directed(points, target) == .02
    # The closest target vertex is much closer than the long segment midpoint;
    # a nearest-midpoint shortcut or a radius U alone would miss the answer.
    assert np.linalg.norm(points[0]-target[3]) < .3
    assert np.linalg.norm(points[0]-.5*(target[0]+target[1])) > 400.


def test_near_touch_and_reversed_full_branches_match_original_distance():
    target = np.array([[0., 0.], [.9, 0.], [.900001, .08], [.900003, -.05], [1., .035]])
    points = np.array([[.900002, .02], [.7, 1e-12], [.900001, .080000000001],
                       [.9000000001, -.000000001]])
    for ordered_points in (points, points[::-1]):
        for ordered_target in (target, target[::-1]):
            assert comparison._directed(ordered_points, ordered_target) == brute_reference(ordered_points, ordered_target)


@pytest.mark.parametrize("scale,translation", [
    (1e-12, np.array([0., 0.])), (1., np.array([0., 0.])),
    (1e6, np.array([-3.2e8, 7.1e8])), (1e-3, np.array([3.2e9, -2.1e9])),
    (1e-140, np.array([0., 0.])),
])
def test_translated_and_scaled_inputs_keep_identical_brute_values(scale, translation):
    rng = np.random.default_rng(52)
    target = scale*rng.normal(size=(93, 2)).cumsum(axis=0)+translation
    points = scale*rng.normal(size=(151, 2))*3+translation
    assert comparison._directed(points, target) == brute_reference(points, target)


def test_tree_queries_are_exact_and_pruning_really_removes_distant_segments(monkeypatch):
    real_tree = comparison.cKDTree
    selected = []
    class Tree:
        def __init__(self, points):
            self.tree = real_tree(points)
        def query(self, points, **kwargs):
            assert kwargs["eps"] == 0.
            return self.tree.query(points, **kwargs)
        def query_ball_point(self, points, radius, **kwargs):
            assert kwargs["eps"] == 0.
            result = self.tree.query_ball_point(points, radius, **kwargs)
            selected.extend(map(len, result))
            return result
    monkeypatch.setattr(comparison, "cKDTree", Tree)
    target = np.column_stack((np.linspace(-10., 10., 201), np.zeros(201)))
    points = np.column_stack((np.linspace(-9., 9., 193), np.full(193, .001)))
    assert comparison._directed(points, target) == brute_reference(points, target)
    assert max(selected) < 5 and len(selected) == len(points)


def test_complete_quadratic_curve_sampling_and_hausdorff_bounds_are_unchanged(monkeypatch):
    first = np.array([[-1., 0., .9, .99, 1.], [0., -.01, -.015, .02, .035]])
    second = first.copy(); second[0, 3] -= .01
    pruned = comparison.curve_distance(first, second, maximum_spacing_m=.002)
    monkeypatch.setattr(comparison, "_directed", brute_reference)
    original = comparison.curve_distance(first, second, maximum_spacing_m=.002)
    assert pruned == original


def test_empty_candidate_query_falls_back_for_every_source_point(monkeypatch):
    real_tree = comparison.cKDTree
    class EmptyCandidates:
        def __init__(self, points):
            self.tree = real_tree(points)
        def query(self, points, **kwargs):
            return self.tree.query(points, **kwargs)
        def query_ball_point(self, points, radius, **kwargs):
            return [[] for _ in points]
    monkeypatch.setattr(comparison, "cKDTree", EmptyCandidates)
    target = np.array([[0., 0.], [1., .2], [.7, .8]])
    points = np.column_stack((np.linspace(-.2, 1.2, 257), np.zeros(257)))
    assert comparison._directed(points, target) == brute_reference(points, target)


@pytest.mark.parametrize("failure", [ValueError, OverflowError, FloatingPointError])
def test_unsuitable_tree_falls_back_to_brute_without_changing_answer(monkeypatch, failure):
    def unusable(*args, **kwargs):
        raise failure("intentional tree arithmetic failure")
    monkeypatch.setattr(comparison, "cKDTree", unusable)
    target = np.array([[0., 0.], [1., .5], [.8, 2.]])
    points = np.array([[.4, .1], [.6, -.3]])
    assert comparison._directed(points, target) == brute_reference(points, target)


def test_nonfinite_or_overflow_inputs_do_not_enter_tree_path(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("invalid tree input was not sent to brute fallback")
    monkeypatch.setattr(comparison, "cKDTree", forbidden)
    cases = [(np.array([[np.inf, 0.]]), np.array([[0., 0.], [1., 0.]])),
             (np.array([[0., 0.]]), np.array([[-1e200, 0.], [1e200, 0.]]))]
    with np.errstate(over="ignore", invalid="ignore"):
        for points, target in cases:
            # This preserves legacy arithmetic for invalid arrays; such arrays
            # are NOT valid native-data acceptance inputs or a PASS fixture.
            assert comparison._directed(points, target) == brute_reference(points, target)


def test_empty_source_keeps_zero_and_target_without_segments_keeps_rejection():
    assert comparison._directed(np.empty((0, 2)), np.array([[0., 0.], [1., 0.]])) == 0.
    with pytest.raises(ValueError):
        comparison._directed(np.array([[0., 0.]]), np.array([[0., 0.]]))
