"""Replay ONE recorded rejected PW2 trial read-only, not a resumed trajectory."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import importlib
import types

import h5py
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'sloshing_visualization/src'))
from sloshing.pinned_wetting import free_boundary as fb
from sloshing.pinned_wetting.free_boundary_remesh import rezone
from sloshing.pinned_wetting.free_boundary_run import identity, REL_OUTPUT
from sloshing.pinned_wetting.free_boundary_verify import (
    _quadratic_curve_intersections, _solid_contact_events, determinant_minima, verifier_dependencies)
from sloshing.pinned_wetting.mission import doctor, job, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path)
    parser.add_argument('--rejected-index', type=int, default=0)
    parser.add_argument('--source-snapshot', action='store_true',
                        help='Execute hash-verified numerical source snapshot; never resume/modify native')
    args = parser.parse_args()
    if args.rejected_index < 0:
        parser.error('rejected-index must be nonnegative')
    cfg = yaml.safe_load((ROOT / 'mission/pinned_wetting/CONFIG.yaml').read_text())
    with job(ROOT, 'PW2 one recorded rejected-trial replay; native read-only', cfg['resources']) as (_, heartbeat):
        before = sha256(args.native)
        with h5py.File(args.native, 'r') as h:
            saved = json.loads(h.attrs['numerical_identity'])
            numerical = fb
            remap = rezone
            if args.source_snapshot:
                snapshot=args.native.parent/'source_snapshot'
                for name,expected in saved['source_hashes'].items():
                    if sha256(snapshot/name)!=expected:
                        raise ValueError('Numerical snapshot hash mismatch: '+name)
                # Isolated import namespace: no replacement of live production
                # modules, and no source copying/patching into the native run.
                alias='pw2_verified_replay_snapshot'
                for package,path in ((alias,snapshot),(alias+'.pinned_wetting',snapshot/'pinned_wetting')):
                    module=types.ModuleType(package);module.__path__=[str(path.resolve())]
                    sys.modules[package]=module
                numerical=importlib.import_module(alias+'.pinned_wetting.free_boundary')
                remap=importlib.import_module(alias+'.pinned_wetting.free_boundary_remesh').rezone
                current=doctor(ROOT)
                locked=json.loads((ROOT/'mission/pinned_wetting/solve_to_animation_v2/CONTRACT.json').read_text())['physical_contract']
                if (saved['physical_contract']!=locked or saved['model_id']!=fb.MODEL_ID
                        or numerical.MODEL_ID!=fb.MODEL_ID or saved.get('synthetic') is not False
                        or any(current['versions'].get(name)!=version for name,version in saved['libraries'].items())):
                    raise ValueError('Snapshot physics/model/libraries are incompatible with this replay.')
            c = numerical.Controls(**saved['controls'])
            if not args.source_snapshot and identity(ROOT, c)[0] != saved:
                raise ValueError('Current numerical sources/settings differ; HEAD alone is not compared.')
            rejection = dict(h['rejected_trials'][f'{args.rejected_index:08d}'].attrs)
            if bool(rejection.get('accepted', True)):
                raise ValueError('Selected journal entry is not explicitly rejected.')
            if 'dt_s' not in rejection:
                raise ValueError('This entry is a remesh failure, not a recorded physical-step trial.')
            n = int(rejection['from_step']); dt = float(rejection['dt_s'])
            group = h['states'][f'{n:08d}']; t = float(group.attrs['time_s'])
            if not bool(group.attrs.get('accepted', False)):
                raise ValueError('Selected starting state is not explicitly accepted.')
            if abs(float(rejection['time_s'])-t)>1e-12 or not 0.<dt<=c.dt+1e-12:
                raise ValueError('Rejected-trial provenance does not match its accepted start.')
            row = json.loads(group.attrs['diagnostics'])
            X, v = group['geometry'][:], group['velocity'][:]
            ids = h['topology/interface_nodes'][:]
            orientation = h['topology/orientation'][:]
            triangles = h['topology/triangles'][:]
            initial_geometry = h['topology/initial_geometry'][:]
        reference = numerical.initial_mesh(c)
        if not (np.array_equal(reference.t, triangles) and np.array_equal(reference.p, initial_geometry)
                and np.array_equal(numerical.interface_nodes(reference), ids)):
            raise ValueError('Regenerated initial topology/geometry differs from native.')
        mesh = replace(reference, doflocs=X)
        native_velocity = v.copy()
        relative_J = np.min(numerical.quadratic_minima(mesh)[0] / numerical.quadratic_minima(reference)[0])
        quality_due = c.rezone_strategy == 'quality_optimized' and relative_J < .01
        due = c.rezone_interval_s > 0 and (
            t >= row.get('last_rezone_time_s', 0.) + c.rezone_interval_s - 1e-11 or quality_due)
        report = {'scope':'ONE_RECORDED_REJECTED_TRIAL_REPLAY_NOT_ACCEPTED', 'accepted':False,
            'native':str(args.native.resolve()), 'native_sha256':before,
            'rejected_index':args.rejected_index, 'from_step':n, 'time_start_s':t, 'dt_s':dt,
            'recorded_rejection':str(rejection['error']), 'rezone_due':bool(due),
            'execution_HEAD':doctor(ROOT)['execution_HEAD'], 'driver_sha256':sha256(__file__),
            'numerical_source_mode':'verified_native_source_snapshot' if args.source_snapshot else 'current_compatible_sources',
            'numerical_source_hashes':saved['source_hashes'],
            'independent_geometry_dependencies':verifier_dependencies(),
            'configurations':[], 'new_native_states_created':0}
        parent = ROOT / REL_OUTPUT / 'rejected_trial_replays'
        parent.mkdir(parents=True, exist_ok=True)
        stem = f'{args.native.parent.name}_{args.rejected_index}_{before[:8]}_'
        output = Path(tempfile.mkdtemp(prefix=stem, dir=parent))
        report['artifact_directory'] = str(output)
        transfer = None
        phase = 'rezoning' if due else 'nonlinear_step'
        try:
            if due:
                mesh, v, transfer = remap(reference, mesh, v, c, orientation)
            phase = 'nonlinear_step'
            next_mesh, next_v, p, diagnostics = numerical.step(mesh, v, dt, c, orientation)
        except (ValueError, RuntimeError) as error:
            report.update(replayed_failure_phase=phase, replayed_solver_rejection=str(error))
        else:
            geometry = {'native_accepted_start':X, 'physical_step_start':mesh.p,
                        'trial_midpoint':.5*(mesh.p+next_mesh.p),
                        'trial_endpoint':next_mesh.p}
            # Preserve the actual candidate BEFORE a separate audit can fail.
            npz = output / 'candidate.npz'
            transferred_arrays = {f'transfer_{key}':value for key,value in (transfer or {}).items()
                                  if isinstance(value,np.ndarray)}
            np.savez_compressed(npz, **geometry, initial_geometry=reference.p, triangles=triangles,
                                interface_ids=ids, native_start_velocity=native_velocity,
                                physical_step_start_velocity=v, end_velocity=next_v,
                                pressure_midpoint=p, dt_s=dt, time_start_s=t, accepted=False,
                                **transferred_arrays)
            report.update(candidate_npz=str(npz), candidate_npz_sha256=sha256(npz),
                          nonlinear_step_completed=True, diagnostics=diagnostics,
                          transfer_metadata={key:value for key,value in (transfer or {}).items()
                                             if not isinstance(value,np.ndarray)})
            try:
                for name, coordinates in geometry.items():
                    candidate = replace(reference, doflocs=coordinates)
                    report['configurations'].append({'name':name,
                        'min_signed_J_m2':float(np.min(determinant_minima(candidate, orientation))),
                        'P2_xy_drift_m':float(np.max(abs(coordinates[:,ids[-1]]-reference.p[:,ids[-1]]))),
                        'quadratic_events':_quadratic_curve_intersections(coordinates, ids),
                        'solid_contact_events':_solid_contact_events(coordinates, ids)})
            except (ValueError, RuntimeError) as error:
                report['independent_geometry_audit_error'] = str(error)
        report['native_unchanged'] = before == sha256(args.native)
        if not report['native_unchanged']:
            raise RuntimeError('Native changed during read-only replay')
        write_json(output / 'report.json', report)
        heartbeat()
        print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
