#!/usr/bin/env python3
"""Check the actual Newton direction against the AD derivative of its merit."""
import argparse,json,runpy,sys
from pathlib import Path
import numpy as np
from scipy.sparse import diags,save_npz
from scipy.sparse.linalg import splu
from geometric_jacobian import pattern,assemble
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);a=p.parse_args();c=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py';sys.argv=[str(runner),'--case',str(c),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];x=ns['state'].copy();n=len(x)//6;norm=solver.getOption('normalizeStates');scale=np.r_[np.full(3*n,norm['U']),np.full(n,norm['p']),np.full(n,norm['T']),np.full(n,norm['nuTilda'])]
def residual(y):solver.solver.updateOFFields(np.ascontiguousarray(y));return solver.getResiduals().copy()
r=residual(x);J=assemble(x,scale,residual,pattern(c,n));D=1/np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-30);save_npz(c/'jacobian.npz',J);np.savez(c/'linearization.npz',state=x,residual=r,scale=scale,row_scale=D)
delta=splu((diags(D)@J).tocsc()).solve(-D*r);g=np.zeros_like(x);solver.setStates(x);solver.solverAD.calcJacTVecProduct('states','stateVar',x,'residuals','residual',np.ascontiguousarray(D*D*r),g)
m0=np.linalg.norm(D*r);slope=float(g@delta/m0**2);fd_gradient=J.T@(D*D*r)
report=dict(raw_l2=float(np.linalg.norm(r)),raw_max=float(abs(r).max()),merit=float(m0),ad_normalized_merit_slope=slope,fd_matrix_normalized_merit_slope=float(fd_gradient@delta/m0**2),weighted_gradient_relative_difference=float(np.linalg.norm(g-fd_gradient)/np.linalg.norm(g)),linear_relative_residual=float(np.linalg.norm(J@delta+r)/np.linalg.norm(r)),trials=[],direction={})
for name,a0,b0 in [('U',0,3*n),('p',3*n,4*n),('T',4*n,5*n),('nuTilda',5*n,6*n)]:
 dx=(scale*delta)[a0:b0];idx=int(np.argmax(abs(dx)));report['direction'][name]=dict(max_abs_change=float(abs(dx).max()),cell=idx//3 if name=='U' else idx,value=float(x[a0+idx]))
for alpha in [1.,.5,.1,.01,.001,1e-4,1e-6,1e-8]:
 y=x+alpha*scale*delta
 if y[3*n:4*n].min()<=1 or y[4*n:5*n].min()<=50 or y[5*n:].min()<0:
  report['trials'].append(dict(alpha=alpha,status='outside_positive_domain'));continue
 rt=residual(y);m=np.linalg.norm(D*rt);report['trials'].append(dict(alpha=alpha,merit=float(m),merit_ratio=float(m/m0),raw_l2=float(np.linalg.norm(rt)),measured_normalized_merit_slope=float((m*m-m0*m0)/(2*alpha*m0*m0))))
solver.setStates(x);np.savez(c/'direction.npz',delta=delta,ad_merit_gradient=g);(c/'direction_report.json').write_text(json.dumps(report,indent=2)+'\n');print('DIRECTION_REPORT '+json.dumps(report),flush=True)
