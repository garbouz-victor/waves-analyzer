#!/usr/bin/env python3
"""Check wall-memory/marker CSVs. Not a CFD or physical-contact validator.

Only Python's standard library is required. Synthetic fixtures belong in tests,
not in scientific output. See OUTPUT_SCHEMA.md for column definitions.
"""
from __future__ import annotations
import argparse
import bisect
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


class CheckError(ValueError):
    """A data-contract violation."""


def number(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise CheckError(f'{label}: expected finite number, got {value!r}') from exc
    if not math.isfinite(result):
        raise CheckError(f'{label}: NaN/inf is forbidden')
    return result


def table(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    try:
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise CheckError(f'{path}: missing or duplicate column names')
            missing = set(required) - set(reader.fieldnames)
            if missing:
                raise CheckError(f'{path}: missing columns {sorted(missing)}')
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise CheckError(f'{path}: {exc}') from exc
    if not rows:
        raise CheckError(f'{path}: empty data')
    if any(None in row or any(row.get(k) is None for k in required) for row in rows):
        raise CheckError(f'{path}: malformed row')
    return rows


def close(actual: float, expected: float, tol: float, label: str) -> None:
    if abs(actual - expected) > tol:
        raise CheckError(f'{label}: {actual:.17g} != {expected:.17g}, tolerance={tol:g}')


def verify(directory: Path, *, z1: float, z2: float, a: float = 1.0,
           t_end: float = 5.0, tolerance_m: float = 1e-9,
           time_tolerance_s: float = 1e-10) -> dict[str, Any]:
    """Validate invariant logic only; caller must separately validate native data."""
    for key, value in {'z1': z1, 'z2': z2, 'a': a, 't_end': t_end,
                       'tolerance_m': tolerance_m, 'time_tolerance_s': time_tolerance_s}.items():
        number(value, key)
    if z2 <= z1 or a <= 0 or t_end <= 0 or tolerance_m <= 0 or time_tolerance_s <= 0:
        raise CheckError('Require z2>z1, a>0, t_end>0 and positive tolerances')
    directory = Path(directory)
    rows = table(directory/'wetting_history.csv',
                 ('time_s','R_L_m','R_R_m','reach_L_m','reach_R_m','H_L_m','H_R_m','state_id'))
    times: list[float] = []
    heights: list[dict[str,float]] = []
    state_ids: list[str] = []
    previous = {'L': z1, 'R': z2}
    for idx, row in enumerate(rows):
        t = number(row['time_s'], f'history[{idx}].time_s')
        if idx == 0:
            close(t, 0.0, time_tolerance_s, 'initial time')
        elif t <= times[-1]:
            raise CheckError(f'history[{idx}]: time is not strictly increasing')
        if t < -time_tolerance_s or t > t_end + time_tolerance_s:
            raise CheckError(f'history[{idx}]: time outside case horizon')
        if not row['state_id'].strip():
            raise CheckError(f'history[{idx}]: empty state_id')
        current = {}
        for side, initial in (('L', z1), ('R', z2)):
            r = number(row[f'R_{side}_m'], f'R_{side}[{idx}]')
            reach = number(row[f'reach_{side}_m'], f'reach_{side}[{idx}]')
            h = number(row[f'H_{side}_m'], f'H_{side}[{idx}]')
            if idx == 0:
                for label, value in (('R',r),('reach',reach),('H',h)):
                    close(value, initial, tolerance_m, f'initial {side} {label}')
            if reach + tolerance_m < r:
                raise CheckError(f'{side} step {idx}: interval reach below endpoint R')
            if h + tolerance_m < previous[side]:
                raise CheckError(f'{side} step {idx}: wet front decreased')
            expected = max(previous[side], reach)
            close(h, expected, tolerance_m, f'{side} step {idx}: running maximum')
            current[side] = h
            previous[side] = h
        times.append(t); heights.append(current); state_ids.append(row['state_id'])
    close(times[-1], t_end, time_tolerance_s, 'final history time')

    def index_at(t: float) -> int:
        k = bisect.bisect_left(times, t)
        choices = [j for j in (k-1,k) if 0 <= j < len(times)]
        if not choices:
            raise CheckError(f'No accepted history near t={t}')
        j = min(choices, key=lambda j: abs(times[j]-t))
        if abs(times[j]-t) > time_tolerance_s:
            raise CheckError(f'Marker creation t={t} not linked to accepted history')
        return j

    cat_rows = table(directory/'marker_catalog.csv',
                     ('marker_id','side','height_m','created_at_s','state_id'))
    catalog: dict[str,dict[str,Any]] = {}
    for row in cat_rows:
        marker = row['marker_id'].strip()
        if not marker or marker in catalog:
            raise CheckError('Empty or reused marker_id')
        side = row['side']
        if side not in ('L','R'):
            raise CheckError(f'{marker}: invalid side')
        z = number(row['height_m'], f'{marker}.height')
        t = number(row['created_at_s'], f'{marker}.created_at')
        j = index_at(t)
        if row['state_id'] != state_ids[j]:
            raise CheckError(f'{marker}: source state_id differs from accepted history')
        close(z, heights[j][side], tolerance_m, f'{marker}: height vs wet maximum at creation')
        catalog[marker] = dict(side=side, height_m=z, created_at_s=t)
    for marker, side, initial in (('P1','L',z1),('P2','R',z2)):
        if marker not in catalog:
            raise CheckError(f'Missing initial marker {marker}')
        item = catalog[marker]
        if item['side'] != side:
            raise CheckError(f'{marker}: wrong wall')
        close(item['height_m'], initial, tolerance_m, f'{marker}: initial height')
        close(item['created_at_s'],0.0,time_tolerance_s,f'{marker}: created at zero')
    for side, initial in (('L',z1),('R',z2)):
        old = initial
        later = sorted((v for k,v in catalog.items() if k not in ('P1','P2') and v['side']==side),
                       key=lambda v:v['created_at_s'])
        for item in later:
            if item['created_at_s'] <= time_tolerance_s or item['height_m'] <= old + tolerance_m:
                raise CheckError(f'{side}: new marker is not a strictly higher later record')
            old = item['height_m']

    track_rows = table(directory/'marker_tracks.csv', ('time_s','marker_id','side','height_m'))
    frames: dict[float,set[str]] = {}
    for row in track_rows:
        t = number(row['time_s'],'track time')
        if t < -time_tolerance_s or t > t_end + time_tolerance_s:
            raise CheckError('Track time outside case horizon')
        marker = row['marker_id']
        if marker not in catalog:
            raise CheckError(f'Track references missing marker {marker}')
        item = catalog[marker]
        if t + time_tolerance_s < item['created_at_s']:
            raise CheckError(f'{marker}: appears before creation')
        if row['side'] != item['side']:
            raise CheckError(f'{marker}: moved to another wall')
        close(number(row['height_m'],f'{marker} track height'), item['height_m'], tolerance_m,
              f'{marker}: immutable position')
        frame = frames.setdefault(t,set())
        if marker in frame:
            raise CheckError(f'{marker}: duplicate in frame {t}')
        frame.add(marker)
    close(min(frames),0.0,time_tolerance_s,'first track frame')
    close(max(frames),t_end,time_tolerance_s,'last track frame')
    for t, actual in frames.items():
        expected = {m for m,v in catalog.items() if v['created_at_s'] <= t+time_tolerance_s}
        if actual != expected:
            raise CheckError(f'Frame {t}: missing persistent markers {sorted(expected-actual)}')
    return {
      'status':'WETTING_CSV_INVARIANTS_PASS',
      'scientific_status':'NOT_VALIDATED_BY_THIS_CHECK',
      'limitations':['Contact heights must be checked against native solver data.',
                     'This does not validate PDE, film connectivity implementation, images or convergence.'],
      'history_rows':len(rows),'track_frames':len(frames),'markers':len(catalog),
      'final_H_L_m':heights[-1]['L'],'final_H_R_m':heights[-1]['R'],
      'x_wall_L_m':-a,'x_wall_R_m':a
    }


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--z1',type=float,required=True)
    p.add_argument('--z2',type=float,required=True)
    p.add_argument('--a',type=float,default=1.0)
    p.add_argument('--t-end',type=float,default=5.0)
    p.add_argument('--tolerance-m',type=float,default=1e-9)
    p.add_argument('--time-tolerance-s',type=float,default=1e-10)
    args=p.parse_args()
    try:
        result=verify(args.directory,z1=args.z1,z2=args.z2,a=args.a,t_end=args.t_end,
                      tolerance_m=args.tolerance_m,time_tolerance_s=args.time_tolerance_s)
    except CheckError as exc:
        print(json.dumps({'status':'FAIL','reason':str(exc)},ensure_ascii=False,indent=2))
        return 1
    print(json.dumps(result,ensure_ascii=False,indent=2)); return 0

if __name__=='__main__':
    sys.exit(main())
