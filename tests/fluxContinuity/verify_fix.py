#!/usr/bin/env python3
"""Check the production library at the independently located flow reversal."""
import argparse,json,runpy,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);a=p.parse_args()
case=a.case.resolve();reference=json.loads(a.reference.read_text());runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py'
sys.argv=[str(runner),'--case',str(case),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart']
ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];state=ns['state'].copy();n=len(state)//6
index=3*n+reference['owner'];root=reference['pressure_shift_to_zero_mach'];rows=[]
for row in reference['samples']:
 if row['mode']!='linear_viscous_only' or row['epsilon_pressure_Pa']>1e-4:continue
 eps=row['epsilon_pressure_Pa'];res=[]
 for sign in [-1,1]:
  x=state.copy();x[index]+=root+sign*eps;solver.setStates(x);res.append(solver.getResiduals().copy())
 jump=float(abs(res[1][4*n:5*n]-res[0][4*n:5*n]).max())
 assert np.isclose(jump,row['energy_residual_jump_max'],rtol=1e-4,atol=1e-9),(jump,row)
 rows.append(dict(epsilon_pressure_Pa=eps,energy_residual_jump_max=jump))
assert rows[-1]['energy_residual_jump_max']<1e-6
assert rows[-1]['energy_residual_jump_max']<rows[0]['energy_residual_jump_max']*1e-3
solver.setStates(state)
report=dict(passed=True,samples=rows,loaded_libraries=ns['loaded'],scope='Production native frozen-state residual continuity; ADR library built and loaded, not a gradient validation.')
(case/'verified_fix.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
