#!/usr/bin/env python3
"""Independent central-FD checks of scaled AD transpose products at a frozen state."""
import argparse,json,runpy,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);a=p.parse_args();case=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py'
sys.argv=[str(runner),'--case',str(case),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];x=ns['state'].copy();n=len(x)//6
normal=solver.getOption('normalizeStates');scales=np.concatenate([np.full(3*n,normal['U']),np.full(n,normal['p']),np.full(n,normal['T']),np.full(n,normal['nuTilda'])])
regions={'U':slice(0,3*n),'p':slice(3*n,4*n),'T':slice(4*n,5*n),'nuTilda':slice(5*n,6*n)}
rng=np.random.default_rng(20261104);seeds={};products={};rows=[]
for name,region in regions.items():
 w=np.zeros(len(x));w[region]=rng.normal(size=len(w[region]));w/=np.linalg.norm(w);seeds[name]=w
 product=np.zeros(len(x));solver.setStates(x);solver.solverAD.calcJacTVecProduct('states','stateVar',x,'residuals','residual',w,product);products[name]=product
for name,region in regions.items():
 v=np.zeros(len(x));v[region]=rng.uniform(-1,1,size=len(v[region]))
 if name=='nuTilda':v[region]*=np.minimum(1,abs(x[region])/scales[region])
 for h in ([1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,3e-9,1e-9,3e-10,1e-10] if name=='T' else [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]):
  solver.setStates(x+h*scales*v);rp=solver.getResiduals().copy();solver.setStates(x-h*scales*v);rm=solver.getResiduals().copy();fd=(rp-rm)/(2*h)
  solver.setStates(x+h*scales*v);ap=np.zeros_like(rp);solver.solverAD.getResiduals(ap)
  solver.setStates(x-h*scales*v);am=np.zeros_like(rm);solver.solverAD.getResiduals(am);fd_ad=(ap-am)/(2*h)
  for output in regions:
   av=float(products[output]@v);fv=float(seeds[output]@fd)
   rows.append(dict(input=name,output=output,h=h,ad=av,fd=fv,relative_error=abs(av-fv)/max(abs(av),abs(fv),1e-12),fd_adr=float(seeds[output]@fd_ad),adr_relative_error=abs(av-float(seeds[output]@fd_ad))/max(abs(av),abs(float(seeds[output]@fd_ad)),1e-12)))
solver.setStates(x);rn=solver.getResiduals().copy();ra=np.zeros_like(rn);solver.solverAD.getResiduals(ra)
checks=[]
for i in regions:
 for o in regions:
  candidates=[r for r in rows if r['input']==i and r['output']==o];best=min(candidates,key=lambda r:r['relative_error']);best=best.copy();best['best_adr_relative_error']=min(r['adr_relative_error'] for r in candidates);checks.append(best)
report=dict(scope='Random block directions; not a full Jacobian or steady gradient verification',normalization='AD product includes normalizeStates; FD perturbs x+h*scales*v',native_ad_residual_relative_difference=float(np.linalg.norm(ra-rn)/np.linalg.norm(rn)),native_ad_residual_max_difference=float(abs(ra-rn).max()),best_by_block=checks,samples=rows,loaded_libraries=ns['loaded'])
(case/'derivative_report.json').write_text(json.dumps(report,indent=2)+'\n');print('DERIVATIVE_SUMMARY '+json.dumps({k:v for k,v in report.items() if k!='samples'}),flush=True)
