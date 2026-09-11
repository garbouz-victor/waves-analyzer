"""Synthetic contract/export fixtures only. No Navier-Stokes computation."""
import copy
import csv
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

MISSION=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(MISSION/'tools'))
from contract_guard import read_json, contract_errors
from check_exports import check_exports, CONTACT_COLUMNS
BASE=read_json(MISSION/'CONTRACT.json')

def save_rows(path,header,data):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=header);w.writeheader();w.writerows(data)

def table(path):
    with Path(path).open(newline='') as f:
        r=csv.DictReader(f);return r.fieldnames,list(r)

def mutate_table(path,edit):
    h,r=table(path);edit(r);save_rows(path,h,r)

def fixture(path):
    """Artificial 11-step export for checking file invariants, NOT physical data."""
    path=Path(path);a=1.;z=math.tan(math.radians(2));depth=10.
    (path/'resolved_case.json').write_text(json.dumps({'physical_contract':BASE['physical_contract'],
        'synthetic_fixture':True}))
    values=[-z,-.01,.01,.005,.015,.012,.008,.004,0.,-.002,-.004]
    contacts=[];intervals=[];curves=[];H=-z
    for i,left in enumerate(values):
        t=i*.5;H=max(H,left)
        contacts.append(dict(zip(CONTACT_COLUMNS,[i,t,left,z,H,z,0.,0.,0.,0.])))
        intervals.extend([{'step':i,'side':'L','z_bottom_m':-depth,'z_top_m':H},
                          {'step':i,'side':'R','z_bottom_m':-depth,'z_top_m':z}])
        for j,x in enumerate([-a,-.5,0.,.5,a]):
            height=left+(x+a)/(2*a)*(z-left)
            curves.append({'time_s':t,'vertex_id':j,'x_m':x,'z_m':height})
    save_rows(path/'contacts.csv',CONTACT_COLUMNS,contacts)
    save_rows(path/'retained_intervals.csv',('step','side','z_bottom_m','z_top_m'),intervals)
    save_rows(path/'interface_samples.csv',('time_s','vertex_id','x_m','z_m'),curves)
    markers=[{'marker_id':'P1','side':'L','x_m':-a,'z_m':-z,'created_at_s':0.,'source_step':0},
             {'marker_id':'P2','side':'R','x_m':a,'z_m':z,'created_at_s':0.,'source_step':0},
             {'marker_id':'P3','side':'L','x_m':-a,'z_m':.01,'created_at_s':1.5,'source_step':2}]
    save_rows(path/'marker_catalog.csv',('marker_id','side','x_m','z_m','created_at_s','source_step'),markers)

class ContractTests(unittest.TestCase):
    def test_same_contract(self):
        self.assertFalse(contract_errors(BASE,{'physical_contract':copy.deepcopy(BASE['physical_contract'])}))
    def test_runtime_lineage_not_physics(self):
        self.assertFalse(contract_errors(BASE,{'physical_contract':BASE['physical_contract'],
            'provenance':{'HEAD':'another'},'numerics':{'mesh':999}}))
    def test_missing_contract(self):self.assertTrue(contract_errors(BASE,{}))
    def test_numeric_int_float_equality(self):
        c=copy.deepcopy(BASE['physical_contract']);c['geometry']['half_width_m']=1
        self.assertFalse(contract_errors(BASE,{'physical_contract':c}))
    def test_boolean_not_numeric(self):
        c=copy.deepcopy(BASE['physical_contract']);c['geometry']['half_width_m']=True
        self.assertTrue(contract_errors(BASE,{'physical_contract':c}))
    def test_duplicate_json_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{"x":1,"x":2}')
            with self.assertRaises(ValueError):read_json(p)
    def test_nan_json_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json';p.write_text('{"x":NaN}')
            with self.assertRaises(ValueError):read_json(p)

def contract_mutation(path,value):
    def test(self):
        c=copy.deepcopy(BASE['physical_contract']);target=c
        for k in path[:-1]:target=target[k]
        target[path[-1]]=value
        self.assertTrue(contract_errors(BASE,{'physical_contract':c}))
    return test
for name,path,value in [
    ('right_slip',('boundaries','right','type'),'Navier'),
    ('left_no_slip',('boundaries','left','type'),'no_slip'),
    ('bottom_slip',('boundaries','bottom','type'),'Navier'),
    ('surface_tension',('material','surface_tension_N_m'),.072),
    ('angle_tuning',('initial','alpha_deg'),.02),
    ('slip_tuning',('boundaries','left','slip_length_m'),.01),
    ('viscosity_tuning',('material','kinematic_viscosity_m2_s'),1.),
    ('linear_relabelled',('target_equations',),'linear_on_z0'),
    ('prewet',('wetting','initial_global_precursor'),True),
    ('unknown_physics_key',('boundaries','right','slip_length_m'),.5),
    ('nan_slip',('boundaries','left','slip_length_m'),float('nan'))]:
    setattr(ContractTests,'test_reject_'+name,contract_mutation(path,value))

class ExportTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.path=Path(self.temp.name);fixture(self.path)
    def tearDown(self):self.temp.cleanup()
    def test_valid_synthetic_export_is_not_cfd_proof(self):
        r=check_exports(self.path,BASE);self.assertTrue(r['passed'],r);self.assertFalse(r['cfd_verified'])
    def test_no_P4_required(self):self.assertTrue(check_exports(self.path,BASE)['passed'])
    def test_nonmonotonic_x_not_rejected_as_linear_slope(self):
        def change(r):
            selected=[v for v in r if float(v['time_s'])==1.]
            selected[1]['x_m']=.2;selected[2]['x_m']=.1
        mutate_table(self.path/'interface_samples.csv',change)
        self.assertTrue(check_exports(self.path,BASE)['passed'])

def export_mutation(filename,edit):
    def test(self):
        mutate_table(self.path/filename,edit)
        r=check_exports(self.path,BASE);self.assertFalse(r['passed'],r)
    return test
for name,filename,edit in [
    ('moving_right','contacts.csv',lambda r:r[2].update(R_right_m=.0)),
    ('coating_shrinks','contacts.csv',lambda r:r[4].update(H_left_m=-.01)),
    ('artificial_H','contacts.csv',lambda r:r[3].update(H_left_m=.9)),
    ('right_H_moved','contacts.csv',lambda r:r[2].update(H_right_m=.04)),
    ('right_slip_velocity','contacts.csv',lambda r:r[3].update(right_w_max_m_s=.01)),
    ('bottom_moves','contacts.csv',lambda r:r[3].update(bottom_speed_max_m_s=.01)),
    ('left_penetrates','contacts.csv',lambda r:r[3].update(left_normal_max_m_s=.01)),
    ('nan','contacts.csv',lambda r:r[3].update(R_left_m='nan')),
    ('missing_final','contacts.csv',lambda r:r.pop()),
    ('duplicate_time','contacts.csv',lambda r:r[3].update(time_s=r[2]['time_s'])),
    ('noncontiguous_steps','contacts.csv',lambda r:r[3].update(step=9)),
    ('dry_hole','retained_intervals.csv',lambda r:r[5].update(z_bottom_m=0)),
    ('missing_retained','retained_intervals.csv',lambda r:r.pop()),
    ('duplicate_retained','retained_intervals.csv',lambda r:r.append(dict(r[0]))),
    ('fake_P2_curve','interface_samples.csv',lambda r:r[14].update(z_m=.0)),
    ('curved_initial','interface_samples.csv',lambda r:r[2].update(z_m=.03)),
    ('missing_snapshot','interface_samples.csv',lambda r:r.__setitem__(slice(None),[v for v in r if float(v['time_s'])!=2.5])),
    ('wrong_left_trace','interface_samples.csv',lambda r:r[10].update(z_m=.02)),
    ('initial_marker_moves','marker_catalog.csv',lambda r:r[1].update(z_m=.0)),
    ('marker_lookahead','marker_catalog.csv',lambda r:r[2].update(created_at_s=.1)),
    ('fake_marker_height','marker_catalog.csv',lambda r:r[2].update(z_m=.011)),
    ('right_record','marker_catalog.csv',lambda r:r[2].update(side='R',x_m=1)),
    ('duplicate_marker','marker_catalog.csv',lambda r:r.append(dict(r[0]))),
]:
    setattr(ExportTests,'test_reject_'+name,export_mutation(filename,edit))

if __name__=='__main__':unittest.main()
