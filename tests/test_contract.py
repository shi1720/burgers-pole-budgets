"""Domain checks needed by the heat-series and contour certificates."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pole_budgets import *

class ContractTests(unittest.TestCase):
 def setUp(self): ctx.prec=160
 def test_time_and_height_domains(self):
    for tau in ('0','-1'):
     for y in (0,1):
      with self.assertRaises(ValueError):support(0,y,tau)
      with self.assertRaises(ValueError):real_point_upper(acb(1),y,tau)
    with self.assertRaises(ValueError):real_point_upper(acb(1),-1,1)
    with self.assertRaises(ValueError):Potential('1').tail(-1)
 def test_truncation_and_mode_types(self):
    for n in (-1,1.5,True):
     with self.assertRaises(ValueError):Potential('1',N=n)
    for k in (.5,True):
     with self.assertRaises(ValueError):FinitePotential({k:1},'1')
 def test_supplied_strip_bound(self):
    p=Potential('0.5');z=acb(0,1)
    with self.assertRaises(ValueError):p.value(z,'.5')
    self.assertTrue(p.value(z,'1').overlaps(p.value(z)))
 def test_contour_resource_contract(self):
    p=Potential('1');rect=('0','2','-1','1')
    for r in [rect[::-1],('0','0','-1','1')]:
     with self.assertRaises(ValueError):radius_bounds(p,r)
    for tol in (-.1,0,1,2,float('nan')):
     with self.assertRaises(ValueError):radius_bounds(p,rect,rtol=tol)
    with self.assertRaises(ValueError):radius_bounds(p,rect,max_boxes=3)
    with self.assertRaises(ValueError):radius_bounds(p,rect,real='False')
    with self.assertRaises(ValueError):certified_count(p,rect,depth=-1)

if __name__=='__main__':unittest.main()
