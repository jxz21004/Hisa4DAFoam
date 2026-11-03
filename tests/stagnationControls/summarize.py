#!/usr/bin/env python3
"""Compare completed single-factor controls using the same residual definitions."""
from pathlib import Path
import json,re,statistics,os
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/hisa-mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3]/'dafoam_case/validation';w=root/'hisa_stagnation_controls';rows=[]
fig,ax=plt.subplots(figsize=(8,4))
for name,c in [('high_order_base',root/'hisa_viscous_fix_L3_2000/L3/base'),('first_order',w/'first_order/L3/base'),('high_order_cfl5',w/'high_order_cfl5/L3/base'),('high_order_from_first',w/'high_order_from_first/L3/base')]:
 if not (c/'primal_report.json').exists():continue
 report=json.loads((c/'primal_report.json').read_text());s=(c/'primal.log').read_text();h=[]
 for raw in re.findall(r'GMRES iteration:\s*0\s+Residual:\s*(.*?)\s+Krylov steps:',s):h.append(max(float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',raw)))
 row=dict(case=name,report=report,counts=json.loads((c/'iteration_counts.json').read_text()) if (c/'iteration_counts.json').exists() else {'krylov_steps_total':sum(int(v) for v in re.findall(r'Krylov steps:\s*(\d+)',s))},outer_samples={str(i):h[i-1] for i in [200,500,1000,1500,2000] if i<=len(h)},last100_mean=statistics.mean(h[-100:]),last100_min=min(h[-100:]),last100_max=max(h[-100:]))
 # High-order cases share the operator and initial state. First-order needs its own initial residual.
 initial=root/'hisa_viscous_fix_L3_2000/initial_probe/final_state.npz' if name!='first_order' else w/'first_order/initial_probe/final_state.npz'
 if initial.exists() and name != "high_order_from_first":
  r0=np.load(initial)['residual'];r=np.load(c/'final_state.npz')['residual'];row['fixed_definition_relative_l2']=float(np.linalg.norm(r)/np.linalg.norm(r0))
 rows.append(row);ax.semilogy(range(1,len(h)+1),h,label=name)
ax.set_xlabel('Outer iteration');ax.set_ylabel('Max HiSA initial normalized residual');ax.grid(True,alpha=.3);ax.legend();fig.tight_layout();fig.savefig(w/'convergence.png',dpi=160)
(w/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
for row in rows:print(json.dumps({k:v for k,v in row.items() if k not in ['report','counts']}),row['report']['raw_residual_max'])
