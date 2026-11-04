#!/usr/bin/env python3
"""Summarize the bounded high-order correction experiments without convergence claims."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/hisa-mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[3];w=root/'dafoam_case/validation/hisa_limited_fix';reports=[json.loads((w/name/'newton_report.json').read_text()) for name in ['geometric','geometric_continuation']]
history=[r for report in reports for r in report['history']];raw=[history[0]['raw_l2_before']]+[r['raw_l2_after'] for r in history if 'raw_l2_after' in r];alphas=[r['alpha'] for r in history if 'raw_l2_after' in r]
fig,axes=plt.subplots(1,2,figsize=(10,4));axes[0].semilogy(range(len(raw)),raw,'o-');axes[0].set_ylabel('Fixed-definition residual L2');axes[1].semilogy(range(1,len(alphas)+1),alphas,'o-');axes[1].set_ylabel('Accepted step length')
for ax in axes:ax.set_xlabel('Accepted correction');ax.grid(True,alpha=.3)
fig.tight_layout();fig.savefig(w/'newton_progress.png',dpi=160)
report=dict(accepted_corrections=len(alphas),raw_l2_initial=raw[0],raw_l2_final=reports[-1]['raw_l2'],raw_max_final=reports[-1]['raw_max'],converged=reports[-1]['accepted'],solver='Experimental coupled finite-difference Newton with direct sparse LU and line search; high-order equations retained',seconds_total=sum(r['seconds'] for r in reports),row_scaling='fixed within each stage, reset at continuation start',stages=reports)
(w/'newton_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='stages'},indent=2))
