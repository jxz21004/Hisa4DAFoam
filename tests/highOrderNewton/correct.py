#!/usr/bin/env python3
"""Experimental sparse quasi-Newton corrections with a true-residual line search."""
import argparse,json,runpy,sys,time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,diags
from scipy.sparse.linalg import splu
from petsc4py import PETSc
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);p.add_argument('--iterations',type=int,default=3);p.add_argument('--geometric',action='store_true');a=p.parse_args();case=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py';sys.argv=[str(runner),'--case',str(case),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];x=ns['state'].copy();n=len(x)//6;normal=solver.getOption('normalizeStates');scale=np.concatenate([np.full(3*n,normal['U']),np.full(n,normal['p']),np.full(n,normal['T']),np.full(n,normal['nuTilda'])]);solver.setOption('adjPartDerivFDStep',{'State':1e-7});solver.updateDAOption();history=[];row_scale=None;start=time.monotonic()
if a.geometric:
 from geometric_jacobian import pattern,assemble
 structure=pattern(case,n)
else:solver.solver.runColoring()
def residual(x):
 solver.solver.updateOFFields(np.ascontiguousarray(x));return solver.getResiduals().copy()
for iteration in range(a.iterations):
 r=residual(x);record=dict(iteration=iteration,raw_l2_before=float(np.linalg.norm(r)),raw_max_before=float(abs(r).max()))
 if abs(r).max()<1e-8:record['status']='converged';history.append(record);break
 if a.geometric:J=assemble(x,scale,residual,structure)
 else:
  mat=PETSc.Mat().create(PETSc.COMM_WORLD);solver.solver.calcdRdWT(0,mat);ptr,col,data=mat.getValuesCSR();J=csr_matrix((data.copy(),col.copy(),ptr.copy()),shape=mat.getSize()).T.tocsc();mat.destroy();J.eliminate_zeros()
 if row_scale is None:row_scale=1/np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-30)
 v=np.random.default_rng(7).normal(size=len(x));v[5*n:]*=np.minimum(1.,x[5*n:]/scale[5*n:]);v/=np.linalg.norm(v);h=1e-6
 fd=(residual(x+h*scale*v)-residual(x-h*scale*v))/(2*h);record['jacobian_directional_relative_error']=float(np.linalg.norm(J@v-fd)/np.linalg.norm(fd));record['jacobian_nnz']=J.nnz
 (case/'matrix_check.json').write_text(json.dumps(record,indent=2)+'\n');print('MATRIX_CHECK '+json.dumps(record),flush=True)
 try:factor=splu((diags(row_scale)@J).tocsc());delta=factor.solve(-row_scale*r)
 except RuntimeError as err:record.update(status='factorization_failed',error=str(err));history.append(record);break
 record['linear_relative_residual']=float(np.linalg.norm(J@delta+r)/np.linalg.norm(r));dx=scale*delta
 active=(x[5*n:]<1e-12)&(dx[5*n:]<0);dx[5*n:][active]=0
 floors=np.r_[np.full(3*n,-np.inf),np.full(n,1.),np.full(n,50.),np.zeros(n)];shrinking=(dx<0)&np.isfinite(floors);alpha=1.
 if np.any(shrinking):alpha=min(alpha,float(np.min(.95*(x[shrinking]-floors[shrinking])/-dx[shrinking])))
 before=np.linalg.norm(row_scale*r);trials=[];accepted=False
 for _ in range(25):
  if alpha<1e-10:break
  trial=x+alpha*dx;rt=residual(trial);merit=np.linalg.norm(row_scale*rt);trials.append(dict(alpha=alpha,merit=float(merit),raw_l2=float(np.linalg.norm(rt))))
  if np.isfinite(merit) and merit<before*(1-1e-4*alpha):x=trial;accepted=True;break
  alpha*=.5
 record.update(status='accepted' if accepted else 'line_search_failed',alpha=alpha,merit_before=float(before),trials=trials)
 if accepted:record.update(raw_l2_after=float(np.linalg.norm(rt)),raw_max_after=float(abs(rt).max()))
 history.append(record);(case/'newton_history.json').write_text(json.dumps(history,indent=2)+'\n');np.savez(case/'newton_state.npz',state=x,residual=residual(x));print('NEWTON_STEP '+json.dumps(record),flush=True)
 if not accepted:break
r=residual(x);report=dict(accepted=bool(abs(r).max()<1e-8),raw_l2=float(np.linalg.norm(r)),raw_max=float(abs(r).max()),seconds=time.monotonic()-start,history=history,scope='Experimental colored-FD quasi-Newton; matrix consistency recorded, fixed row-scaled merit, positive-state safeguards; no steady gradient claim.')
(case/'newton_report.json').write_text(json.dumps(report,indent=2)+'\n');np.savez(case/'newton_state.npz',state=x,residual=r);print('NEWTON_REPORT '+json.dumps(report),flush=True)
