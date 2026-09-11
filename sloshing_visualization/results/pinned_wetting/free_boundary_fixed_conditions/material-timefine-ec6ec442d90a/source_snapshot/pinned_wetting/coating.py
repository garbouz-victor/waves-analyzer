"""Zero-volume, lower-connected retained liquid, independent of any renderer."""
from dataclasses import asdict, dataclass, field
import math


@dataclass(frozen=True)
class Marker:
    marker_id: str
    side: str
    height_m: float
    created_at_s: float
    state_id: str


@dataclass
class Coating:
    run_id: str
    bottom: float
    H: list
    markers: list = field(default_factory=list)
    accepted_step: int = 0
    time: float = 0.0

    @classmethod
    def initial(cls, run_id, bottom, z1, z2):
        if not all(map(math.isfinite, [bottom, z1, z2])) or not bottom < z1 < z2:
            raise ValueError("Invalid initial geometry")
        sid = f"{run_id}:0"
        return cls(run_id, bottom, [z1, z2],
                   [Marker("P1", "L", z1, 0., sid), Marker("P2", "R", z2, 0., sid)])

    def accept(self, step, time_s, R, reach=None):
        """Call only AFTER the solver accepts the complete interval.

        A solver with internal extrema must supply qualified reach, not a
        rejected stage value. Endpoint sampling is explicit in the reference.
        """
        reach = list(R if reach is None else reach)
        if (step != self.accepted_step + 1 or time_s <= self.time or
                len(R) != 2 or len(reach) != 2 or
                not all(map(math.isfinite, [time_s, *R, *reach])) or
                any(p < r or r < self.bottom for p, r in zip(reach, R))):
            raise ValueError("Invalid accepted contact interval")
        self.H = [max(h, p) for h, p in zip(self.H, reach)]
        self.accepted_step, self.time = step, time_s

    def liquid_present(self, side, z):
        """The entire interval is wet; no mesh or plotting threshold can make holes."""
        if side not in ("L", "R") or not math.isfinite(z):
            raise ValueError("Invalid wall query")
        return self.bottom <= z <= self.H[0 if side == "L" else 1]

    def record_peak(self, marker_id, side, height, peak_time, peak_state_id):
        if marker_id in {m.marker_id for m in self.markers} or side not in ("L", "R"):
            raise ValueError("Marker ID reused or side invalid")
        prior = max(m.height_m for m in self.markers if m.side == side)
        if (not all(map(math.isfinite, [height, peak_time])) or
                height <= prior or height > self.H[0 if side == "L" else 1] or
                not 0 <= peak_time <= self.time or
                not peak_state_id.startswith(self.run_id + ":")):
            raise ValueError("Unqualified record")
        # Creation time is confirmation time, not a look-ahead to the peak.
        self.markers.append(Marker(marker_id, side, height, self.time, peak_state_id))

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        value["markers"] = [Marker(**m) for m in value["markers"]]
        result = cls(**value)
        if (len(result.H) != 2 or not all(map(math.isfinite, result.H)) or
                min(result.H) < result.bottom or
                len({m.marker_id for m in result.markers}) != len(result.markers)):
            raise ValueError("Corrupt retained state")
        return result
