import sys,json,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pole_budgets import *
ROOT=Path(__file__).resolve().parents[1];ctx.prec=160
rng=np.random.default_rng(1720406)
class Perturbed:
 def __init__(self,p,c): self.p=p;self.c=c
 def value(self,z,Y=None):
  return self.p.value(z,Y)+sum((acb(v.real,v.imag)*(-self.p.T*k*k).exp()*(acb(0,k)*z).exp() for k,v in self.c.items()),acb(0))
rows=[]
for tau,Y in [('0.125','1'),('0.5','1'),('1','2')]:
 p=Potential(tau);rect=('0','2','-'+Y,Y);base=certified_count(p,rect)['count']
 for real in (False,True):
  budget=radius_bounds(p,rect,real=real)['lower_float']*.45
  for j in range(8):
   c={k:complex(*rng.normal(size=2)) for k in range(-4,5)}
   if real:
    c[0]=complex(c[0].real,0)
    for k in range(1,5):c[-k]=c[k].conjugate()
   norm=sum(map(abs,c.values()));c={k:v*budget/norm for k,v in c.items()}
   r=certified_count(Perturbed(p,c),rect)
   assert r['status']=='certified' and r['count']==base
   rows.append({'tau':tau,'Y':Y,'real':real,'draw':j,'budget':budget,'count':r['count']})
# Explicit two-mode necessity at tau=Y=1: E0 and E1 determine the hull.
t=arb(1);aa=(-t).exp()*t.cosh();bb=(-t).exp()*t.sinh()
th=((1-aa*aa)/(bb*bb)).sqrt().atan();w=acb((1+aa*aa)/2,bb*bb*th.tan()/2)
dual=(w.real*th.cos()+w.imag*th.sin())/support(th,1,1)
single=((w.real/aa)**2+(w.imag/bb)**2).sqrt()
geometry={'target':str(w),'theta':str(th),'dual':str(dual),'two_mode_cost':'1 (exact algebraic construction)','best_single_pair':str(single),'verified_upper':str(real_point_upper(w,1,1))}
# Invariance under representation and precision: N 24/40/56 and 96/160/224 bits.
repeats=[]
for prec in (96,160,224):
 ctx.prec=prec
 for N in (24,40,56):
  p=Potential('0.5',N=N);r=certified_count(p,('0','2','-1','1'))
  assert r['status']=='certified' and r['count']==2
  q=radius_bounds(p,('0','2','-1','1'),real=True)
  assert q['converged'];repeats.append({'bits':prec,'N':N,'count':r['count'],'radius':q})
(ROOT/'results/validation.json').write_text(json.dumps({'seed':1720406,'perturbation_checks':rows,'two_mode':geometry,'precision_truncation':repeats},indent=2))
print('48 perturbation count checks, 9 precision/truncation runs; two-mode single cost',single)
