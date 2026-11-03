#!/usr/bin/env python3
"""Frozen-state continuity probe. Does not advance any primal time steps."""
import argparse,os,sys,runpy,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);a=p.parse_args()
case=a.case.resolve();lib=a.library_dir.resolve();dump=case/'faces.csv'
os.environ['HISA_PROBE_CSV']=str(dump)
runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py'
sys.argv=[str(runner),'--case',str(case),'--library-dir',str(lib),'--probe','--restart']
ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];state=ns['state'].copy();n=len(state)//6
saved=np.load(case/'final_state.npz')['residual'].copy()
def evaluate(x):
 solver.setStates(x);r=solver.getResiduals().copy();f=np.genfromtxt(dump,delimiter=',',names=True)
 return r,f
r0,f0=evaluate(state)
np.savez(case/'baseline_probe.npz',state=state,residual=r0)
jump=np.sqrt(sum((f0['UL'+d]-f0['UR'+d])**2 for d in 'xyz'))
candidates=np.flatnonzero(jump>1.)
candidates=candidates[np.argsort(abs(f0['mach'][candidates]))]
selected=None
for face in candidates[:30]:
 cell=int(f0['owner'][face]);idx=3*n+cell
 x=state.copy();x[idx]+=1.;_,f=evaluate(x)
 dm=f['mach'][face]-f0['mach'][face]
 if abs(dm)<1e-14:continue
 estimate=-f0['mach'][face]/dm
 if abs(estimate)>0.01*state[idx]:continue
 lo=min(0.,estimate*2);hi=max(0.,estimate*2)
 def sample(dp):
  x=state.copy();x[idx]+=dp;return evaluate(x)
 _,fl=sample(lo);_,fh=sample(hi)
 if fl['mach'][face]*fh['mach'][face]<0:
  selected=(int(face),cell,idx,lo,hi,fl['mach'][face]);break
assert selected is not None,'No nearby pressure-induced flow reversal found'
face,cell,idx,lo,hi,ml=selected
for _ in range(45):
 mid=.5*(lo+hi);_,fm=sample(mid);mm=fm['mach'][face]
 if mm*ml>0:lo=mid;ml=mm
 else:hi=mid
root=.5*(lo+hi)
results=[]
for mode in ['upwind','linear_viscous_only']:
 if mode=='upwind':os.environ.pop('HISA_PROBE_LINEAR_VISCOUS',None)
 else:os.environ['HISA_PROBE_LINEAR_VISCOUS']='1'
 for eps in [1e-2,1e-4,1e-6,1e-8]:
  rm,fm=sample(root-eps);rp,fp=sample(root+eps)
  diff=rp-rm
  results.append(dict(mode=mode,epsilon_pressure_Pa=eps,mach_minus=float(fm['mach'][face]),mach_plus=float(fp['mach'][face]),Up_jump=float(np.linalg.norm([fp['Up'+d][face]-fm['Up'+d][face] for d in 'xyz'])),mass_flux_jump=float(fp['massFlux'][face]-fm['massFlux'][face]),convective_energy_flux_jump=float(fp['energyFlux'][face]-fm['energyFlux'][face]),residual_jump_max=float(abs(diff).max()),energy_residual_jump_max=float(abs(diff[4*n:5*n]).max()),energy_residual_owner_jump=float(diff[4*n+cell])))
os.environ.pop('HISA_PROBE_LINEAR_VISCOUS',None)
solver.setStates(state)
report=dict(face=face,owner=cell,neighbour=int(f0['neighbour'][face]),baseline_mach=float(f0['mach'][face]),baseline_velocity_reconstruction_jump=float(jump[face]),baseline_owner_pressure=float(state[idx]),pressure_shift_to_zero_mach=root,repeat_residual_max_difference=float(abs(r0-saved).max()),samples=results,scope='Frozen L3 state; only one cell pressure varied. Diagnostic linear viscous interpolation is not a production fix or a convergence claim.')
(case/'continuity_report.json').write_text(json.dumps(report,indent=2)+'\n')
print('CONTINUITY_REPORT '+json.dumps(report),flush=True)
