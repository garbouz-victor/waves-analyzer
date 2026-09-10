"""Synthetic data only. Tests never simulate or certify a fluid trajectory."""
from pathlib import Path
import csv
import importlib.util
import tempfile
import unittest

MODULE=Path(__file__).resolve().parents[1]/'tools/check_wetting_csv.py'
spec=importlib.util.spec_from_file_location('pw_csv_checker',MODULE)
checker=importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)

H=('time_s','R_L_m','R_R_m','reach_L_m','reach_R_m','H_L_m','H_R_m','state_id')
C=('marker_id','side','height_m','created_at_s','state_id')
T=('time_s','marker_id','side','height_m')

class WettingTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.history=[
          [0,-.1,.1,-.1,.1,-.1,.1,'s0'],
          [1,.2,-.05,.2,.05,.2,.1,'s1'],
          [2,0,.08,.1,.08,.2,.1,'s2']]
        self.catalog=[['P1','L',-.1,0,'s0'],['P2','R',.1,0,'s0'],['P3','L',.2,1,'s1']]
        self.tracks=[]
        for t in (0,1,2):
            for marker,side,z,created,_ in self.catalog:
                if t>=created: self.tracks.append([t,marker,side,z])
    def tearDown(self): self.tmp.cleanup()
    def runcheck(self):
        for name,cols,rows in [('wetting_history.csv',H,self.history),
                               ('marker_catalog.csv',C,self.catalog),
                               ('marker_tracks.csv',T,self.tracks)]:
            with (self.root/name).open('w',newline='') as h:
                w=csv.writer(h); w.writerow(cols); w.writerows(rows)
        return checker.verify(self.root,z1=-.1,z2=.1,t_end=2)
    def test_valid(self): self.assertEqual(self.runcheck()['status'],'WETTING_CSV_INVARIANTS_PASS')
    def test_not_physical_certification(self):
        self.assertEqual(self.runcheck()['scientific_status'],'NOT_VALIDATED_BY_THIS_CHECK')
    def test_decrease(self):
        self.history[2][5]=.1
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_invented_rise(self):
        self.history[2][5]=.4
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_marker_moves(self):
        self.tracks[-1][3]=.3
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_marker_disappears(self):
        self.tracks=self.tracks[:-1]
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_nan(self):
        self.history[1][1]='nan'
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_short_horizon(self):
        self.history=self.history[:2]
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_initial_point_wrong(self):
        self.catalog[1][2]=.12
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_duplicate_marker(self):
        self.catalog.append(self.catalog[-1][:])
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_no_new_record_is_valid(self):
        self.history=[[0,-.1,.1,-.1,.1,-.1,.1,'s0'],[2,-.2,0,-.2,0,-.1,.1,'s2']]
        self.catalog=self.catalog[:2]
        self.tracks=[r for r in self.tracks if r[1] in ('P1','P2')]
        self.assertEqual(self.runcheck()['markers'],2)
    def test_right_new_record_does_not_move_p2(self):
        self.history[2][2]=.3; self.history[2][4]=.3; self.history[2][6]=.3
        self.catalog.append(['P4','R',.3,2,'s2']); self.tracks.append([2,'P4','R',.3])
        self.assertEqual(self.runcheck()['final_H_R_m'],.3)
    def test_internal_peak_retained(self):
        self.history[1][1]=.1  # endpoint below peak but qualified reach=.2
        self.assertEqual(self.runcheck()['final_H_L_m'],.2)
    def test_reach_below_R(self):
        self.history[1][3]=.1
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_same_time_twice(self):
        self.history[2][0]=1
        with self.assertRaises(checker.CheckError): self.runcheck()
    def test_source_id_mismatch(self):
        self.catalog[-1][-1]='made-up'
        with self.assertRaises(checker.CheckError): self.runcheck()

if __name__=='__main__': unittest.main()
