#!/usr/bin/env python3
"""Gate use of the corrected operator on primal equivalence and AD parity."""
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('fixed',type=Path);p.add_argument('original',type=Path);a=p.parse_args()
f=np.load(a.fixed/'parity_fields.npz');b=np.load(a.original/'parity_fields.npz');rn=f['residual_native'];ra=f['residual_adr']
unchanged=float(abs(rn-b['residual_native']).max());parity=float(np.linalg.norm(rn-ra)/np.linalg.norm(rn));d=json.loads((a.fixed/'derivative_report.json').read_text());worst=max(r['relative_error'] for r in d['best_by_block'])
report=dict(primal_max_difference=unchanged,native_adr_relative_difference=parity,worst_best_directional_relative_error=worst,passed=unchanged<1e-8 and parity<1e-12 and worst<1e-4)
(a.fixed/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));assert report['passed'],report
