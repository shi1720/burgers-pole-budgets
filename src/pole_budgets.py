"""Validated contour counts and coefficient budgets for periodic heat potentials.
All decisions use Arb balls. Floating point only selects test directions/segments.
"""
from dataclasses import dataclass
from fractions import Fraction
import math
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, acb, ctx


def A(x):
    if isinstance(x, Fraction): return arb(str(x.numerator))/x.denominator
    return arb(x)


def maxball(xs):
    # Outward interval containing the maximum of a finite family.
    return arb(max(x.lower() for x in xs)).union(arb(max(x.upper() for x in xs)))


def amplification(y,tau):
    """Exact infinite-mode operator norm exp(max_k(k*y-tau*k^2))."""
    y=A(y); tau=A(tau)
    if not tau>0: raise ValueError('tau must be positive')
    if not y>=0: raise ValueError('y must be nonnegative')
    j=max(0,int(float((y/(2*tau)).upper())))
    # floor from float is only a proposal; verify the omitted tails decrease.
    K=j+2
    while not 2*tau*K>y: K+=1
    return maxball([(k*y-tau*k*k).exp() for k in range(K+1)])


def support(theta,y,tau):
    """Support of the real-symmetric coefficient l1 unit ball after heat flow."""
    y=A(y); tau=A(tau); theta=A(theta)
    if not tau>0: raise ValueError('tau must be positive')
    if y==0: return abs(theta.cos())
    if not y>0: raise ValueError('strict positive y needed')
    c,s=theta.cos(),theta.sin()
    # h_1(theta) >= exp(-tau)*sinh(y); k beyond K decreases geometrically.
    floor=(-tau).exp()*y.sinh()
    K=max(2,int(float((y/tau).upper()))+2)
    while not (tau*K>y and (-tau*(K+1)**2+(K+1)*y).exp()<floor): K+=1
    vals=[abs(c)]
    vals += [(-tau*k*k).exp()*((k*y).cosh()**2*c*c+(k*y).sinh()**2*s*s).sqrt() for k in range(1,K+1)]
    return maxball(vals)


def choose_theta(w,y,tau):
    """Heuristic proposal; no floating-point result is trusted as a bound."""
    if y<1e-14: return 0.0 if w.real>=0 else math.pi
    K=max(10,int(y/tau)+15)
    k=np.arange(K+1); E=np.exp(-tau*k*k)
    aa=E*np.cosh(k*y); bb=E*np.sinh(k*y)
    def obj(th):
        h=np.max(np.hypot(aa*np.cos(th),bb*np.sin(th)))
        return -(w.real*np.cos(th)+w.imag*np.sin(th))/h
    grid=np.linspace(-math.pi,math.pi,129)
    j=min(range(len(grid)),key=lambda j:obj(grid[j]))
    r=minimize_scalar(obj,bounds=(grid[max(0,j-1)],grid[min(128,j+1)]),method='bounded',options={'xatol':1e-13})
    return float(r.x)


def real_point_upper(w,y,tau):
    """Constructive bound using one conjugate Fourier pair (or constant)."""
    y=A(y); tau=A(tau)
    if not tau>0 or not y>=0: raise ValueError('tau positive and y nonnegative required')
    if y==0:
        if w.imag==0: return abs(w.real)
        return arb.pos_inf()
    vals=[]
    # Any selected mode gives a rigorous feasible witness; no cutoff assumption.
    K=max(8,int(float((y/tau).upper()))+8)
    for k in range(1,K+1):
        vals.append((tau*k*k).exp()*((w.real/(k*y).cosh())**2+(w.imag/(k*y).sinh())**2).sqrt())
    if w.imag==0: vals.append(abs(w.real))
    # Near a nonsmooth support maximum, two modes can be cheaper than one.
    th=A(choose_theta(complex(w),float(y),float(tau)))
    c,s=th.cos(),th.sin(); points=[]
    for k in range(K+1):
        aa=(-tau*k*k).exp()*(k*y).cosh()
        bb=(-tau*k*k).exp()*(k*y).sinh()
        h=((aa*c)**2+(bb*s)**2).sqrt()
        if h>0: points.append((aa*aa*c/h,bb*bb*s/h))
    for i,(px,py) in enumerate(points):
        for qx,qy in points[i+1:]:
            det=px*qy-py*qx
            if not abs(det)>0: continue
            lam=(w.real*qy-w.imag*qx)/det
            mu=(px*w.imag-py*w.real)/det
            if lam>=0 and mu>=0: vals.append(lam+mu)
    return min(v.upper() for v in vals)

@dataclass
class Potential:
    tau: str
    a: str='2'
    N: int=40
    def __post_init__(self):
        self.T=A(self.tau); self.alpha=A(self.a)
        if not self.T>0 or not self.alpha>0 or type(self.N) is not int or self.N<0: raise ValueError('invalid potential')
        self.c=[self.alpha.bessel_i(k)*(-self.T*k*k).exp() for k in range(self.N+1)]
    def tail(self,Y):
        Y=A(Y); n=self.N+1; a=self.alpha; T=self.T
        if not Y>=0: raise ValueError('strip height must be nonnegative')
        q=a/(2*(n+1))*(Y-T*(2*n+1)).exp()
        if not q<1: raise ValueError('increase N for tail bound')
        term=2*(a*a/(4*(n+1))).exp()*(a/2)**n/arb.fac_ui(n)*(-T*n*n+n*Y).exp()
        return term/(1-q)
    def value(self,z,Y=None):
        bound=abs(z.imag).upper()
        Y=bound if Y is None else A(Y).upper()
        if not Y>=bound: raise ValueError('strip height must cover the evaluation box')
        m=acb(z.real.mid(),z.imag.mid())
        r2=z.real.rad()**2+z.imag.rad()**2
        v=acb(self.c[0]); d=acb(0); second=arb(0)
        centered=float(r2)>1e-36
        for k in range(1,self.N+1):
            v+=2*self.c[k]*(k*(m if centered else z)).cos()
            if centered:
                d-=2*k*self.c[k]*(k*m).sin()
                second+=k*k*self.c[k]*(k*Y).cosh()
        e=self.tail(Y).upper()
        if centered:
            v+=d*(z-m)
            e+=(second*r2).upper()
        return v+acb(arb(0,e),arb(0,e))
    def approx(self,z):
        return sum(float(c)*(1 if k==0 else 2*np.cos(k*z)) for k,c in enumerate(self.c))


def box(p,q):
    return acb((A(p[0])*arb.pi()).union(A(q[0])*arb.pi()),A(p[1]).union(A(q[1])))


def certified_count(pot,rect,depth=24):
    """CCW rectangle (horizontal coordinates in multiples of pi); image boxes certify phase increments."""
    x0,x1,y0,y1=map(Fraction,rect)
    if not x0<x1 or not y0<y1: raise ValueError('bad rectangle')
    if type(depth) is not int or depth<0: raise ValueError('depth must be a nonnegative integer')
    pts=[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
    todo=[(pts[j],pts[(j+1)%4],0) for j in reversed(range(4))]
    accepted=[]; calls=0
    while todo:
        p,q,d=todo.pop(); v=pot.value(box(p,q)); calls+=1
        if abs(v)>0:
            accepted.append((p,q)); continue
        if d>=depth: return {'status':'unresolved','evaluations':calls}
        m=((p[0]+q[0])/2,(p[1]+q[1])/2)
        todo.extend([(m,q,d+1),(p,m,d+1)])
    total=arb(0)
    for p,q in accepted:
        f=pot.value(acb(A(p[0])*arb.pi(),A(p[1]))); g=pot.value(acb(A(q[0])*arb.pi(),A(q[1])))
        total+=(g/f).arg()
    turns=total/(2*arb.pi()); n=round(float(turns.mid()))
    if not (turns>arb(n)-arb('0.25') and turns<arb(n)+arb('0.25')):
        return {'status':'unresolved','reason':'phase sum','turns':str(turns)}
    return {'status':'certified','count':n,'segments':len(accepted),'evaluations':calls,'turns':str(turns)}


def radius_bounds(pot,rect,real=False,rtol=.02,max_boxes=20000):
    """Branch and bound for boundary-zero distance, with outward bounds.
    Real lower bounds are dual support bounds; upper bounds use feasible pairs.
    """
    x0,x1,y0,y1=map(Fraction,rect)
    if not x0<x1 or not y0<y1: raise ValueError('bad rectangle')
    if type(real) is not bool: raise ValueError('real must be a boolean')
    if not math.isfinite(rtol) or not 0<rtol<1: raise ValueError('relative tolerance must be between zero and one')
    if type(max_boxes) is not int or max_boxes<4: raise ValueError('at least four boxes are required')
    pts=[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
    upper=arb.pos_inf(); records=[]; evaluations=0
    def evaluate(p,q):
        nonlocal upper,evaluations
        evaluations+=1
        v=pot.value(box(p,q)); mid=((p[0]+q[0])/2,(p[1]+q[1])/2)
        vm=pot.value(acb(A(mid[0])*arb.pi(),A(mid[1]))); Y=max(abs(p[1]),abs(q[1])); ym=abs(mid[1])
        lo=(abs(v).lower()/amplification(Y,pot.T)).lower()
        if real:
            if ym==0:
                up=real_point_upper(-vm,ym,pot.T)
            else:
                th=A(choose_theta(-complex(vm),float(ym),float(pot.T)))
                num=(-v.real*th.cos()-v.imag*th.sin()).lower()
                # h increases in |y|; finite mode cap is checked in support().
                if num>0: lo=max(lo,(num/support(th,Y,pot.T)).lower())
                up=real_point_upper(-vm,ym,pot.T)
        else:
            up=(abs(vm).upper()/amplification(ym,pot.T)).upper()
        if up<upper: upper=up
        return (lo,p,q)
    for j in range(4): records.append(evaluate(pts[j],pts[(j+1)%4]))
    while evaluations<max_boxes:
        j=min(range(len(records)),key=lambda j:records[j][0]); low,p,q=records[j]
        if low>0 and low>=upper*A(str(1-rtol)): break
        m=((p[0]+q[0])/2,(p[1]+q[1])/2)
        if m==p or m==q: break
        records[j]=evaluate(p,m); records.append(evaluate(m,q))
    lower=min(r[0] for r in records)
    return {'lower':str(lower),'upper':str(upper),'lower_float':float(lower),'upper_float':float(upper),
            'evaluations':evaluations,'converged':bool(lower>0 and lower>=upper*A(str(1-rtol)))}

class FinitePotential:
    """Exact finite initial Fourier coefficients indexed by integer mode."""
    def __init__(self,coeff,tau):
        self.T=A(tau)
        if not self.T>0: raise ValueError('tau must be positive')
        if any(type(k) is not int for k in coeff): raise ValueError('Fourier modes must be integers')
        self.c={k:acb(v)*(-self.T*k*k).exp() for k,v in coeff.items()}
    def value(self,z,Y=None):
        return sum((v*(acb(0,k)*z).exp() for k,v in self.c.items()),acb(0))
