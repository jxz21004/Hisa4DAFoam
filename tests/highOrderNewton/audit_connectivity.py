#!/usr/bin/env python3
"""Check assembled columns against identical uncolored forward differences."""
import argparse,runpy,sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from petsc4py import PETSc
from geometric_jacobian import labels
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--library-dir',type=Path,required=True);p.add_argument('--dafoam-library-dir',type=Path);a=p.parse_args();a.dafoam_library_dir=a.dafoam_library_dir.resolve() if a.dafoam_library_dir else None;c=a.case.resolve();runner=Path(__file__).resolve().parents[1]/'wingAero/run_case.py';sys.argv=[str(runner),'--case',str(c),'--library-dir',str(a.library_dir.resolve()),'--probe','--restart'];ns=runpy.run_path(str(runner),run_name='__main__');solver=ns['solver'];loaded_da=sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines() if '/libDASolver' in line});
if a.dafoam_library_dir:
 assert len(loaded_da)==2 and all(Path(v).parent==a.dafoam_library_dir.resolve() for v in loaded_da),loaded_da
x=ns['state'].copy();n=len(x)//6;h=1e-8;normal=solver.getOption('normalizeStates');solver.setOption('adjPartDerivFDStep',{'State':h});solver.updateDAOption();solver.solver.runColoring();m=PETSc.Mat().create(PETSc.COMM_WORLD);solver.solver.calcdRdWT(0,m);ptr,col,data=m.getValuesCSR();J=csr_matrix((data.copy(),col.copy(),ptr.copy()),shape=m.getSize()).T.tocsc();m.destroy()
owner=labels(c/'constant/polyMesh/owner');nei=labels(c/'constant/polyMesh/neighbour');adj=[set() for _ in range(n)]
for u,v in zip(owner,nei):adj[u].add(int(v));adj[v].add(int(u))
def res(y):solver.solver.updateOFFields(np.ascontiguousarray(y));return solver.getResiduals().copy()
r0=res(x);rows=[]
for cell in [9,123,128,100,1567]:
 column=3*cell;y=x.copy();y[column]+=h*normal['U'];fd=(res(y)-r0)/h;assembled=J[:,column].toarray().ravel();missing=np.flatnonzero((abs(fd[5*n:])>1e-6)&(assembled[5*n:]==0))
 examples=[]
 for rowcell in sorted(missing,key=lambda i:abs(fd[5*n+i]),reverse=True)[:5]:
  dist=0 if rowcell==cell else (1 if rowcell in adj[cell] else (2 if any(rowcell in adj[j] for j in adj[cell]) else '>2'))
  examples.append(dict(row_cell=int(rowcell),graph_distance=dist,fd=float(fd[5*n+rowcell]),assembled=float(assembled[5*n+rowcell])))
 rows.append(dict(column_cell=cell,component='Ux',relative_column_error=float(np.linalg.norm(assembled-fd)/max(np.linalg.norm(fd),1e-30)),sa_missing_count=len(missing),examples=examples))
res(x);report=dict(loaded_dafoam={v:hashlib.sha256(Path(v).read_bytes()).hexdigest() for v in loaded_da},fd_step=h,scope='Same forward difference and state scaling; nonzero finite-difference SA coupling absent in assembled column',columns=rows);(c/'connectivity_report.json').write_text(json.dumps(report,indent=2)+'\n');print('CONNECTIVITY_REPORT '+json.dumps(report),flush=True)
