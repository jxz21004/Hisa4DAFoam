#!/usr/bin/env python3
"""Compare diagonal damping of the same frozen high-order Newton system."""
import argparse,runpy,sys,json,time
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz,diags
from scipy.sparse.linalg import splu
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);a=p.parse_args();c=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py';sys.argv=[str(runner),'--case',str(c),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];z=np.load(c/'linearization.npz');x=z['state'];r=z['residual'];scale=z['scale'];D=z['row_scale'];n=len(x)//6;J=load_npz(c/'jacobian.npz');A=diags(D)@J;g=np.load(c/'direction.npz')['ad_merit_gradient'];m0=np.linalg.norm(D*r);history=[];best=None;start=time.monotonic()
def res(y):solver.solver.updateOFFields(np.ascontiguousarray(y));return solver.getResiduals().copy()
assert np.allclose(res(x),r,rtol=1e-12,atol=1e-8)
for damping in [1e-3,1e-2,1e-1]:
 B=A+diags(damping*np.r_[-np.ones(5*n),np.ones(n)]);delta=splu(B.tocsc()).solve(-D*r);dx=scale*delta;sl=float(g@delta/m0**2);alpha=1.;floors=np.r_[np.full(3*n,-np.inf),np.full(n,1.),np.full(n,50.),np.zeros(n)];shrinking=(dx<0)&np.isfinite(floors)
 if np.any(shrinking):alpha=min(alpha,float(np.min(.95*(x[shrinking]-floors[shrinking])/-dx[shrinking])))
 record=dict(damping=damping,ad_merit_slope=sl,trials=[],accepted=False)
 if sl<0:
  for _ in range(20):
   if alpha<1e-10:break
   y=x+alpha*dx;rt=res(y);m=np.linalg.norm(D*rt);record['trials'].append(dict(alpha=alpha,merit_ratio=float(m/m0),raw_l2=float(np.linalg.norm(rt))))
   if np.isfinite(m) and .5*m*m<=.5*m0*m0+1e-4*alpha*sl*m0*m0:
    record.update(accepted=True,alpha=alpha,merit_ratio=float(m/m0),raw_l2=float(np.linalg.norm(rt)),raw_max=float(abs(rt).max()))
    if best is None or m<best[0]:best=(m,y.copy(),rt.copy(),damping)
    break
   alpha*=.5
 history.append(record);print('DAMPING_TRIAL '+json.dumps(record),flush=True);(c/'regularization_trials.json').write_text(json.dumps(history,indent=2)+'\n')
report=dict(initial_raw_l2=float(np.linalg.norm(r)),initial_merit=float(m0),seconds=time.monotonic()-start,trials=history)
if best is not None:np.savez(c/'damped_state.npz',state=best[1],residual=best[2]);report.update(selected_damping=best[3],final_raw_l2=float(np.linalg.norm(best[2])),final_raw_max=float(abs(best[2]).max()))
solver.setStates(x);(c/'regularization_report.json').write_text(json.dumps(report,indent=2)+'\n');print('REGULARIZATION_REPORT '+json.dumps(report),flush=True)
