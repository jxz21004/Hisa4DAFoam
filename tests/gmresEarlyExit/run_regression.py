#!/usr/bin/env python3
"""Three one-step L3 smoke tests of real HiSA GMRES (requires DAFoam environment)."""
import argparse,json,os,re,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--template',type=Path,required=True)
p.add_argument('--library-dir',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
script=Path(__file__).resolve().parents[1]/'wingAero/run_case.py'
env=os.environ.copy(); env['LD_LIBRARY_PATH']=str(a.library_dir.resolve())+':'+env.get('LD_LIBRARY_PATH','')
env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
results=[]
for name,dirs,cycles,rtol in [('early',40,60,0.1),('budget',3,2,1e-30),('single',1,3,1e-30)]:
 c=a.output/name;c.mkdir()
 for d in ['0','constant','system']:shutil.copytree(a.template/d,c/d)
 f=c/'system/controlDict';s=f.read_text()
 for key,value in [('endTime','1'),('writeInterval','1')]:s=re.sub(r'\b'+key+r'\s+[^;]+;',key+' '+value+';',s)
 f.write_text(s)
 f=c/'system/fvSolution';s=f.read_text()
 s=re.sub(r'\bnKrylov\s+\d+;',f'nKrylov {dirs};',s)
 s=re.sub(r'\bmaxIter\s+\d+;',f'maxIter {cycles};',s)
 s=re.sub(r'\bsolverTolRel\s+[^;]+;',f'solverTolRel {rtol} ({rtol} {rtol} {rtol}) {rtol};\n        solverTol 1e-30 (1e-30 1e-30 1e-30) 1e-30;',s)
 f.write_text(s)
 with (c/'primal.log').open('w') as log:
  subprocess.run([sys.executable,str(script),'--case',str(c.resolve()),'--library-dir',str(a.library_dir.resolve())],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 text=(c/'primal.log').read_text()
 pattern=r'GMRES iteration:\s*(\d+)\s+Residual:\s*(.*?)\s+Krylov steps:\s*(\d+)'
 matches=re.findall(pattern,text)
 assert len(matches)==2,(name,matches)
 values=lambda x:[float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?',x)]
 initial=values(matches[0][1]);final=values(matches[-1][1]);steps=int(matches[-1][2]);count=int(matches[-1][0])
 ratios=[f/i for f,i in zip(final,initial) if i>1e-25]
 if name=='early':assert 0<steps<dirs and max(ratios)<rtol,(steps,ratios)
 else:assert steps==dirs*cycles and count==cycles,(name,steps,count)
 results.append(dict(case=name,subspace=dirs,cycles=count,krylov_steps=steps,initial=initial,final=final,max_relative_residual=max(ratios),passed=True))
(a.output/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
