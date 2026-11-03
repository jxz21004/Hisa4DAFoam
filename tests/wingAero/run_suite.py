#!/usr/bin/env python3
"""Run prepared L3/L2 cases and gate steady gradient checks on primal convergence."""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('work',type=Path)
parser.add_argument('--levels',nargs='+',choices=['L3','L2'],default=['L3','L2'])
parser.add_argument('--summary',default='summary.json')
args=parser.parse_args()
work=args.work.resolve(strict=True); libs=work/'build/lib'
script=Path(__file__).with_name('run_case.py').resolve()
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(libs)+':'+env.get('LD_LIBRARY_PATH','')
env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'


def run(case,alpha=2.,restart=False,adjoint=False):
    command=[sys.executable,str(script),'--case',str(case),'--library-dir',str(libs),'--alpha',str(alpha)]
    if restart:command.append('--restart')
    if adjoint:command.extend(['--adjoint','--probe','--restart'])
    log_path=case/('adjoint.log' if adjoint else 'primal.log')
    with log_path.open('w') as log:
        code=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    text=log_path.read_text()
    cycles=[int(n) for n in re.findall(r'GMRES iteration:\s*(\d+)',text)]
    solution=(case/'system/fvSolution').read_text()
    subspace=int(re.search(r'\bnKrylov\s+(\d+)\s*;',solution)[1])
    cfl=[float(value) for value in re.findall(r'^Pseudo Courant No:\s*([^\s]+)',text,re.M)]
    actual_steps=[int(n) for n in re.findall(r'Krylov steps:\s*(\d+)',text)]
    stats={'exit_code':code,'gmres_restart_cycles_total':sum(cycles),'gmres_subspace':subspace,
           'gmres_krylov_steps_total':sum(actual_steps) if actual_steps else subspace*sum(cycles),'gmres_max_restart_cycles_per_solve':max(cycles,default=0),
           'configured_cfl_max':float(re.search(r'\bpseudoCoNumMax\s+([^;]+);',solution)[1]),
           'actual_cfl_max':max(cfl,default=None),'actual_cfl_final':cfl[-1] if cfl else None,
           'temperature_clipping_messages':len(re.findall(r'bound.*[Tt]emperature|[Tt]emperature.*bound|bounding T',text))}
    if not adjoint:(case/'iteration_counts.json').write_text(json.dumps(stats,indent=2)+'\n')
    return code


def level_run(level):
    base=work/level/'base';run(base)
    report_path=base/'primal_report.json'
    if not report_path.exists():
        return {'level':level,'accepted':False,'reason':'Primal process failed before final report','gradient_verified':False}
    baseline=json.loads(report_path.read_text())
    summary={'level':level,'base':baseline,'iterations':json.loads((base/'iteration_counts.json').read_text()),'gradient_verified':False}
    if not baseline['accepted']:
        summary['reason']='Base flow did not pass raw residual gate; no steady-gradient accuracy claim'
        return summary
    times=[p for p in base.iterdir() if p.is_dir() and re.fullmatch(r'\d+(\.\d+)?',p.name)]
    latest=max(times,key=lambda p:float(p.name))
    perturbations={}
    for name,alpha in [('minus',1.999),('plus',2.001)]:
        target=work/level/name;target.mkdir()
        for folder in ('constant','system'):shutil.copytree(base/folder,target/folder)
        shutil.copytree(latest,target/'0')
        if (target/'0/uniform').exists():shutil.rmtree(target/'0/uniform')
        for field in ('pseudoCoField','pseudoCoField.gz'):
            if (target/'0'/field).exists():(target/'0'/field).unlink()
        run(target,alpha,restart=True)
        p=target/'primal_report.json'
        perturbations[name]=json.loads(p.read_text()) if p.exists() else {'accepted':False}
    summary['perturbations']=perturbations
    if not all(p['accepted'] for p in perturbations.values()):
        summary['reason']='Perturbed primal failed; FD steady gradient is not verified'
        return summary
    run(base,adjoint=True)
    p=base/'adjoint_report.json'
    if not p.exists():
        summary['reason']='Adjoint process failed';return summary
    adj=json.loads(p.read_text());summary['adjoint']=adj
    passed=True
    for name,row in adj['functions'].items():
        fd=(perturbations['plus']['functions'][name]-perturbations['minus']['functions'][name])/.002
        error=abs(fd-row['d_dalpha']);row.update(fd=fd,absolute_error=error,relative_error=error/max(abs(fd),1e-30))
        passed &= row['linear_fail']==0 and row['true_linear_relative_residual']<1e-8 and error<=1e-6+1e-3*abs(fd)
    summary['gradient_verified']=bool(passed)
    return summary

results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    futures={pool.submit(level_run,level):level for level in args.levels}
    for future in concurrent.futures.as_completed(futures):
        level=futures[future]
        try:result=future.result()
        except Exception as error:result={'level':level,'error':str(error),'gradient_verified':False}
        results.append(result)
        (work/args.summary).write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(result),flush=True)
