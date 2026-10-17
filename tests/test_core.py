import unittest,sys,math,json
from pathlib import Path
import numpy as np
from scipy.special import iv
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pole_budgets import *

class Tests(unittest.TestCase):
 def setUp(self): ctx.prec=160
 def test_infinite_amplification(self):
    for tau in ('0.03125','0.125','0.5','2'):
     for y in ('0','0.25','1','3'):
      v=amplification(y,tau)
      j=round(float(y)/(2*float(tau)))
      ref=(A(j)*A(y)-A(tau)*j*j).exp()
      self.assertTrue(v.overlaps(ref))
      self.assertTrue(v.contains(max(np.exp(k*float(y)-float(tau)*k*k) for k in range(200))) or abs(float(v)-float(ref))<1e-14)
 def test_exact_counts(self):
    # Heat potential 1 + b*cos(m*z). Zeros are x=(2j+1)pi/m,
    # |y|=acosh(exp(tau*m*m)/b)/m, directly from cos(m*z).
    for m in (1,2,3):
     for tau in ('0.125','0.5'):
      b=.5; y=math.acosh(math.exp(float(tau)*m*m)/b)/m
      p=FinitePotential({0:1,m:'.25',-m:'.25'},tau)
      for fac in (.8,1.2):
       yy=str(y*fac);r=certified_count(p,('0','2','-'+yy,yy))
       self.assertEqual(r['status'],'certified'); self.assertEqual(r['count'],0 if fac<1 else 2*m)
 def test_boundary_refusal(self):
    # 1-exp(-tau)*exp(iz) has a zero at z=i*(-tau).
    p=FinitePotential({0:1,1:-1},'0.5')
    r=certified_count(p,('0','2','-0.5','0.5'),depth=8)
    self.assertEqual(r['status'],'unresolved')
 def test_tail_and_precision(self):
    for tau in ('0.125','0.5'):
     for Y in ('0.5','2'):
      small=Potential(tau,N=12);big=Potential(tau,N=60)
      for x in (.1,1.7,3.1,5.9):
       z=acb(x,Y); self.assertTrue(small.value(z).overlaps(big.value(z)))
      self.assertGreater(float(small.tail(Y)),0)
 def test_real_support_vs_complex(self):
    for tau in ('0.125','0.5','1'):
     for y in ('0.25','1','2'):
      amp=amplification(y,tau)
      for th in (0,.3,1,1.57,2.7): self.assertTrue(support(th,y,tau).lower()<=amp.upper())
 def test_primal_dual(self):
    for tau in ('0.125','0.5','1'):
     for y in ('0.25','1','2'):
      for w in (acb(1,0),acb(0,1),acb(1,1),acb(1,3)):
       theta=A(choose_theta(complex(w),float(y),float(tau)))
       lo=(w.real*theta.cos()+w.imag*theta.sin())/support(theta,y,tau)
       up=real_point_upper(w,y,tau)
       self.assertTrue(up>=lo.lower())
       self.assertLess(float(up/lo),1.00001)
 def test_invalid(self):
    with self.assertRaises(ValueError): Potential('0')
    with self.assertRaises(ValueError): amplification(-1,1)
    with self.assertRaises(ValueError): certified_count(Potential('1'),('2','0','-1','1'))
 def test_radius_symmetry_gap(self):
    p=Potential('0.5');rect=('0','2','-1','1')
    c=radius_bounds(p,rect,max_boxes=3000);r=radius_bounds(p,rect,real=True,max_boxes=3000)
    self.assertTrue(c['converged'] and r['converged'])
    self.assertGreater(r['lower_float'],c['upper_float']*1.6)

if __name__=='__main__': unittest.main(verbosity=2)
