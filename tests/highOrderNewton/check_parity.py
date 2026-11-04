#!/usr/bin/env python3
"""Compare native/ADR fields and residuals at identical state, without iteration."""
import argparse,runpy,sys,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);a=p.parse_args();c=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py';sys.argv=[str(runner),'--case',str(c),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];x=ns['state'].copy();n=len(x)//6;out={};arrays={}
for repeat in range(3):
 solver.setStates(x);rn=solver.getResiduals();ra=np.zeros_like(rn);solver.solverAD.getResiduals(ra)
 out[str(repeat)]={'residual_max_difference':float(abs(rn-ra).max()),'fields':{}}
 for name,typ in [('U','vector'),('p','scalar'),('T','scalar'),('nuTilda','scalar'),('rho','scalar'),('rhoE','scalar'),('nut','scalar'),('alphat','scalar'),('thermo:mu','scalar'),('thermo:alpha','scalar'),('eZero','scalar')]:
  fn=np.zeros(n*(3 if typ=='vector' else 1));fa=fn.copy();solver.solver.getOFField(name,typ,fn);solver.solverAD.getOFField(name,typ,fa)
  out[str(repeat)]['fields'][name]={'max_difference':float(abs(fn-fa).max()),'relative_l2':float(np.linalg.norm(fn-fa)/max(np.linalg.norm(fn),1e-30)),'native_range':[float(fn.min()),float(fn.max())],'adr_range':[float(fa.min()),float(fa.max())]};arrays[name+'_native']=fn;arrays[name+'_adr']=fa
np.savez(c/'parity_fields.npz',**arrays,residual_native=rn,residual_adr=ra)
(c/'parity_report.json').write_text(json.dumps(out,indent=2)+'\n');print('PARITY_REPORT '+json.dumps(out),flush=True)

points=np.zeros(solver.solver.getNLocalPoints()*3);solver.solver.getOFMeshPoints(points)
solver.setVolCoords(points);solver.setStates(x)
rn_new=solver.getResiduals().copy();ra_new=np.zeros_like(rn_new);solver.solverAD.getResiduals(ra_new)
move=dict(native_change=float(np.linalg.norm(rn_new-rn)),adr_change=float(np.linalg.norm(ra_new-ra)),native_adr_max_difference=float(abs(rn_new-ra_new).max()),native_adr_relative_difference=float(np.linalg.norm(rn_new-ra_new)/np.linalg.norm(rn_new)))
(c/'geometry_refresh_report.json').write_text(json.dumps(move,indent=2)+'\n');print('GEOMETRY_REFRESH '+json.dumps(move),flush=True)
