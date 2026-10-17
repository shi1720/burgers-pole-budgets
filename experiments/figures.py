import sys,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from scipy.special import iv
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pole_budgets import *
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
rows=json.loads((ROOT/'results/grid.json').read_text())
# Figure 1: visual location is illustrative; counts come from the validated code.
x=np.linspace(0,2*np.pi,601);y=np.linspace(-2,2,401);z=x[None,:]+1j*y[:,None]
v=sum(iv(k,2)*np.exp(-.5*k*k)*(1 if k==0 else 2*np.cos(k*z)) for k in range(30))
hsv=np.zeros(v.shape+(3,));hsv[:,:,0]=(np.angle(v)/(2*np.pi))%1;hsv[:,:,1]=.7;hsv[:,:,2]=.55+.4*(1-np.exp(-np.abs(v)))
fig,ax=plt.subplots(figsize=(6.5,3.4));ax.imshow(hsv_to_rgb(hsv),origin='lower',extent=[0,2,-2,2],aspect='auto');ax.axhline(1,color='white',ls='--');ax.axhline(-1,color='white',ls='--');ax.set(xlabel=r'$\mathrm{Re}\,z/\pi$',ylabel=r'$\mathrm{Im}\,z$',title=r'Heat potential $\Phi$, $a=2$, $\tau=1/2$: phase and magnitude');ax.text(.04,1.65,'solid frame: count 2\ndashed contour: count 2',color='white');fig.savefig(out/'landscape.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,3.2));colors=['#23628c','#c65b32','#477e48','#7759a5']
for tau,col in zip(['0.125','0.25','0.5','1'],colors):
 r=[r for r in rows if r['tau']==tau];xx=[float(s['Y']) for s in r]
 for ax,key in zip(axs,['complex','real']):
  lo=np.array([s[key]['lower_float'] for s in r]);hi=np.array([s[key]['upper_float'] for s in r]);ax.plot(xx,(lo+hi)/2,'o-',color=col,label=tau);ax.fill_between(xx,lo,hi,color=col,alpha=.25)
for ax,title in zip(axs,['Complex coefficient errors','Errors preserving real data']):ax.set(xlabel='Contour half-height Y',ylabel='Boundary-zero budget',title=title,yscale='log');ax.grid(alpha=.2)
axs[0].legend(title=r'$\tau$',fontsize=8,ncol=2);fig.tight_layout();fig.savefig(out/'budgets.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(6.5,3.3));th=np.linspace(0,2*np.pi,600);tau=1;y=1
for k,col in [(1,'#23628c'),(2,'#8ca8b3'),(3,'#c9d7dc')]:
 a=np.exp(-tau*k*k)*np.cosh(k*y);b=np.exp(-tau*k*k)*np.sinh(k*y);ax.plot(a*np.cos(th),b*np.sin(th),color=col,label=f'mode pair {k}')
a=np.exp(-1)*np.cosh(1);b=np.exp(-1)*np.sinh(1);theta=np.arctan(np.sqrt(1-a*a)/b);p=np.array([a*a,b*b*np.tan(theta)]);w=(p+np.array([1,0]))/2
ax.plot([-1,1],[0,0],color='#c65b32',label='constant mode');ax.plot([1,p[0]],[0,p[1]],'k--',label='exposed hull segment');ax.scatter(*w,c='black');ax.annotate('two-mode target w',w,xytext=(.18,.31),arrowprops={'arrowstyle':'->'});ax.set(xlabel='Real perturbation at z',ylabel='Imaginary perturbation at z',title='Real-data value set: convex hull of Fourier ellipses',aspect='equal');ax.legend(fontsize=8,loc='lower left',ncol=2);fig.savefig(out/'geometry.pdf');plt.close(fig)
# Noise amplification is independent of the nominal solution.
fig,axs=plt.subplots(1,2,figsize=(7,3.0));ys=np.linspace(0,2,400)
for tau,col in zip([.125,.25,.5,1],colors):
 k=np.arange(40);amp=np.exp(np.max(k[:,None]*ys[None,:]-tau*k[:,None]**2,axis=0));j=np.argmax(k[:,None]*ys[None,:]-tau*k[:,None]**2,axis=0)
 axs[0].semilogy(ys,amp,color=col,label=str(tau));axs[1].step(ys,j,color=col,where='mid')
axs[0].set(xlabel='Imaginary height |y|',ylabel=r'$A_\tau(y)$',title='Heat smoothing versus complex amplification');axs[0].legend(title=r'$\tau$',fontsize=8,ncol=2)
axs[1].set(xlabel='Imaginary height |y|',ylabel='Most damaging mode |k|',title='A finite, moving critical frequency');fig.tight_layout();fig.savefig(out/'amplification.pdf');plt.close(fig)
# Automatically generated, outward-rounded result table.
def interval(q):
 lo=q['lower_float'];hi=q['upper_float'];e=int(np.floor(np.log10(hi)));scale=10**(e-3)
 return f'[{np.floor(lo/scale)*scale:.4g}, {np.ceil(hi/scale)*scale:.4g}]'
lines=[]
for r in rows:
 lines.append(f"{r['tau']} & {r['Y']} & {r['count']['count']} & {interval(r['complex'])} & {interval(r['real'])} \\\\")
(ROOT/'paper/results_table.tex').write_text('\n'.join(lines))
print('4 vector figures and result table written')
