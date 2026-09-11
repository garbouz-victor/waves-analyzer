#!/usr/bin/env python3
"""Check exported geometry/contact invariants, NOT hydrodynamic correctness."""
from __future__ import annotations
import argparse
import bisect
import csv
import json
import math
from pathlib import Path
import sys
from contract_guard import DEFAULT_CONTRACT, contract_errors, read_json

CONTACT_COLUMNS = ('step','time_s','R_left_m','R_right_m','H_left_m','H_right_m',
                   'right_u_max_m_s','right_w_max_m_s','bottom_speed_max_m_s','left_normal_max_m_s')
TOL = 1e-12

def rows(path, required):
    with Path(path).open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f'{path.name}: missing/duplicate header columns')
        result = list(reader)
        if not result or any(None in r for r in result):
            raise ValueError(f'{path.name}: empty table or extra unnamed CSV columns')
        return result

def num(row, key):
    v = float(row[key])
    if not math.isfinite(v):
        raise ValueError(f'{key}: nonfinite numeric value')
    return v

def integer(row, key):
    s = row[key]
    if not s or not s.isdecimal():
        raise ValueError(f'{key}: expected nonnegative integer text')
    return int(s)

def check_exports(output: Path, contract: dict) -> dict:
    errors = []
    def require(condition, reason):
        if not condition and reason not in errors:
            errors.append(reason)
    try:
        output = Path(output)
        errors.extend(contract_errors(contract, read_json(output/'resolved_case.json')))
        c = contract['physical_contract']; target = contract['deliverable']
        a, depth = c['geometry']['half_width_m'], c['geometry']['depth_m']
        k = math.tan(math.radians(c['initial']['alpha_deg']))
        z1, z2 = -a*k, a*k
        contact_rows = rows(output/'contacts.csv', CONTACT_COLUMNS)
        contacts = []
        H = [z1,z2]
        prev_time = None
        for index, row in enumerate(contact_rows):
            step = integer(row,'step'); t = num(row,'time_s')
            R = [num(row,'R_left_m'),num(row,'R_right_m')]
            reported_H = [num(row,'H_left_m'),num(row,'H_right_m')]
            reach = [num(row,key) if row.get(key) not in (None,'') else R[i]
                     for i,key in enumerate(('reach_left_m','reach_right_m'))]
            require(step == index, 'noncontiguous accepted steps')
            require(prev_time is None or t > prev_time, 'non-increasing contact time')
            require(target['t_start_s']-TOL <= t <= target['t_end_s']+TOL, 'contact time outside target')
            require(all(-depth <= r for r in R), 'contact below bottom')
            require(all(p >= r-TOL for p,r in zip(reach,R)), 'reach below endpoint')
            if any(p > r+TOL for p,r in zip(reach,R)):
                require(bool(row.get('reach_evidence','').strip()), 'internal reach missing native evidence reference')
            if index == 0:
                require(abs(t-target['t_start_s'])<=TOL, 'missing t0')
                require(abs(R[0]-z1)<=TOL and abs(R[1]-z2)<=TOL, 'wrong initial P1/P2')
                require(all(abs(p-r)<=TOL for p,r in zip(reach,R)), 'initial reach already prewet')
            H = [max(h,p) for h,p in zip(H,reach)]
            require(all(abs(h-v)<=TOL for h,v in zip(H,reported_H)), 'H is not accepted-contact running maximum')
            require(abs(R[1]-z2)<=TOL, 'right current material contact moved')
            require(abs(reported_H[1]-z2)<=TOL, 'right H differs from P2')
            for key in CONTACT_COLUMNS[6:]:
                value = num(row,key)
                require(0 <= value <= TOL, f'boundary trace failure: {key}')
            contacts.append({'step':step,'t':t,'R':R,'H':reported_H})
            prev_time = t
        require(abs(contacts[-1]['t']-target['t_end_s'])<=TOL, 'missing 5s contact endpoint')
        times = [v['t'] for v in contacts]
        # Exact accepted-time matching within output roundoff; not arbitrary temporal interpolation.
        def contact_at(t):
            j = bisect.bisect_left(times,t)
            candidates = [q for q in (j-1,j) if 0<=q<len(times)]
            if not candidates:
                return None
            q = min(candidates,key=lambda x:abs(times[x]-t))
            return contacts[q] if abs(times[q]-t)<=1e-10 else None
        intervals = rows(output/'retained_intervals.csv',('step','side','z_bottom_m','z_top_m'))
        seen = set()
        for row in intervals:
            step,side=integer(row,'step'),row['side']
            require(side in ('L','R') and step<len(contacts), 'invalid retained interval key')
            key = (step,side)
            require(key not in seen, 'duplicate retained interval')
            seen.add(key)
            if side in ('L','R') and step<len(contacts):
                require(abs(num(row,'z_bottom_m')+depth)<=TOL and
                        abs(num(row,'z_top_m')-contacts[step]['H'][side=='R'])<=TOL,
                        'retained interval is not full lower-connected W')
        require(seen=={(i,s) for i in range(len(contacts)) for s in ('L','R')}, 'missing retained intervals')
        interface = rows(output/'interface_samples.csv',('time_s','vertex_id','x_m','z_m'))
        curves = {}
        for row in interface:
            t = num(row,'time_s')
            curves.setdefault(t,[]).append((integer(row,'vertex_id'),num(row,'x_m'),num(row,'z_m')))
        for t in target['snapshot_times_s']:
            require(any(abs(s-t)<=TOL for s in curves), f'missing interface snapshot {t}')
        for t, data in curves.items():
            data.sort()
            require(len(data)>=2 and [p[0] for p in data]==list(range(len(data))), 'curve vertex IDs invalid')
            state = contact_at(t)
            require(state is not None,'interface time lacks accepted contact state')
            require(abs(data[0][1]+a)<=TOL and abs(data[-1][1]-a)<=TOL,'wrong wall positions on curve')
            require(abs(data[-1][2]-z2)<=TOL,'curve does not end at fixed P2')
            if state:
                require(abs(data[0][2]-state['R'][0])<=TOL,'curve left endpoint differs from contact')
            require(all(-a-TOL<=x<=a+TOL and z>=-depth-TOL for _,x,z in data),'curve outside vessel')
            require(all(math.hypot(p[1]-q[1],p[2]-q[2])>0 for p,q in zip(data,data[1:])), 'zero length surface segment')
            if abs(t)<=TOL:
                require(all(abs(z-k*x)<=1e-10 for _,x,z in data),'initial interface is not specified straight line')
        markers=rows(output/'marker_catalog.csv',('marker_id','side','x_m','z_m','created_at_s','source_step'))
        marker_ids = set()
        for row in markers:
            name,side = row['marker_id'],row['side']
            x,z,born = num(row,'x_m'),num(row,'z_m'),num(row,'created_at_s')
            step=integer(row,'source_step')
            require(name not in marker_ids,'duplicate marker ID'); marker_ids.add(name)
            require(side in ('L','R') and step<len(contacts),'invalid marker side/source')
            require(0<=born<=target['t_end_s']+TOL,'marker creation time outside target')
            require(abs(x-(-a if side=='L' else a))<=TOL,'marker not on its wall')
            if name in ('P1','P2'):
                require(side==('L' if name=='P1' else 'R') and abs(z-(z1 if name=='P1' else z2))<=TOL
                        and born==0 and step==0,'initial marker moved/relabelled')
            else:
                require(side=='L','new right marker at fixed P2 forbidden')
                if step<len(contacts):
                    require(born+TOL>=contacts[step]['t'],'marker appears before source event')
                    require(abs(z-contacts[step]['R'][0])<=TOL,'marker height not from accepted contact')
                    require(abs(z-max(v['R'][0] for v in contacts[:step+1]))<=TOL,'marker not record reached by source time')
        require({'P1','P2'}<=marker_ids,'missing initial markers')
    except (OSError,ValueError,KeyError,TypeError,IndexError,OverflowError) as exc:
        errors.append('schema/read error: '+str(exc))
    return {'scope':'EXPORTED_INVARIANTS_ONLY_NOT_PDE','passed':not errors,
            'cfd_verified':False,'errors':errors}

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--contract',type=Path,default=DEFAULT_CONTRACT)
    args=parser.parse_args(argv)
    try:
        result=check_exports(args.output,read_json(args.contract))
    except (ValueError,OSError) as exc:
        result={'scope':'EXPORTED_INVARIANTS_ONLY_NOT_PDE','passed':False,'cfd_verified':False,'errors':[str(exc)]}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0 if result['passed'] else 2

if __name__=='__main__':
    sys.exit(main())
