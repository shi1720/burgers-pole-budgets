import sys,json,time,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pole_budgets import *
import flint
ROOT=Path(__file__).resolve().parents[1]
ctx.prec=160
protocol=json.loads((ROOT/'experiments/protocol.json').read_text())
rows=[]
for tau in protocol['tau']:
 for Y in protocol['Y']:
    p=Potential(tau,a=protocol['a'],N=protocol['N'])
    rect=('0','2','-'+Y,Y); start=time.perf_counter()
    r={'tau':tau,'Y':Y,'N':p.N,'count':certified_count(p,rect)}
    for real in (False,True):
        r['real' if real else 'complex']=radius_bounds(p,rect,real=real,rtol=protocol['rtol'],max_boxes=protocol['max_boxes'])
    r['seconds']=time.perf_counter()-start
    rows.append(r)
    (ROOT/'results/grid.json').write_text(json.dumps(rows,indent=2))
    print(tau,Y,r['count'].get('count'),[(key,round(r[key]['lower_float'],8),r[key]['converged']) for key in ('complex','real')],flush=True)
(ROOT/'results/environment.json').write_text(json.dumps({'python':sys.version,'flint':flint.__version__,'platform':platform.platform(),'precision_bits':ctx.prec},indent=2))
