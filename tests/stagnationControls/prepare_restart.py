#!/usr/bin/env python3
"""Restore high order from the completed first-order solution; reset CFL history."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parents[3];w=root/'dafoam_case/validation/hisa_stagnation_controls';c=w/'high_order_from_first/L3/base';c.mkdir(parents=True,exist_ok=False)
for d in ['constant','system']:shutil.copytree(root/'dafoam_case/validation/hisa_viscous_fix_L3_2000/L3/base'/d,c/d)
shutil.copytree(w/'first_order/L3/base/2000',c/'0')
if (c/'0/uniform').exists():shutil.rmtree(c/'0/uniform')
for name in ['pseudoCoField','pseudoCoField.gz']:
 if (c/'0'/name).exists():(c/'0'/name).unlink()
(w/'high_order_from_first/build').symlink_to(root/'dafoam_case/validation/hisa_viscous_work_fix/build',target_is_directory=True)
print('Run wingAero/run_case.py with --restart, without primalInitCondition.')
