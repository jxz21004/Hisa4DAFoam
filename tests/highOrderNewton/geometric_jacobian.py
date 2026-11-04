"""Independent two-ring cell-block coloring for serial FV residual diagnostics."""
import gzip,re
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,eye

def labels(path):
 path=Path(path)
 text=path.read_text() if path.exists() else gzip.open(str(path)+'.gz','rt').read()
 match=re.search(r'\n\s*(\d+)\s*\(\s*([\d\s]+)\)',text)
 if match is None:raise ValueError('Expected ASCII labelList: '+str(path))
 vals=np.fromstring(match[2],dtype=np.int64,sep=' ');assert len(vals)==int(match[1]);return vals

def pattern(case,n):
 mesh=Path(case)/'constant/polyMesh';owner=labels(mesh/'owner');neighbour=labels(mesh/'neighbour');own=owner[:len(neighbour)]
 A=coo_matrix((np.ones(2*len(own)),(np.r_[own,neighbour],np.r_[neighbour,own])),shape=(n,n)).tocsr()+eye(n,format='csr')
 B=(A@A).astype(bool).astype(float).tocsr();G=(B.T@B).tocsr();G.setdiag(0);G.eliminate_zeros();colors=np.full(n,-1,dtype=int)
 for cell in np.argsort(-np.diff(G.indptr),kind='stable'):
  forbidden=set(colors[G.indices[G.indptr[cell]:G.indptr[cell+1]]]);color=0
  while color in forbidden:color+=1
  colors[cell]=color
 table=np.column_stack([3*np.arange(n)+d for d in range(3)]+[offset*n+np.arange(n) for offset in [3,4,5]])
 return B,colors,table

def assemble(x,scale,rfun,structure,h=1e-8):
 B,colors,table=structure;rows=[];cols=[];data=[]
 for var in range(6):
  for color in range(int(colors.max())+1):
   cells=np.flatnonzero(colors==color);columns=table[cells,var];xp=x.copy();xm=x.copy();xp[columns]+=h*scale[columns];xm[columns]-=h*scale[columns]
   diff=(rfun(xp)-rfun(xm))/(2*h)
   for cell,column in zip(cells,columns):
    adjacent=B.indices[B.indptr[cell]:B.indptr[cell+1]];rr=table[adjacent].ravel();rows.append(rr);cols.append(np.full(len(rr),column));data.append(diff[rr])
  print('GEOMETRIC_JACOBIAN field',var,'colors',int(colors.max())+1,flush=True)
 J=coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(x),len(x))).tocsc();J.eliminate_zeros();return J
