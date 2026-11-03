#!/usr/bin/env python3
"""Compile/link the changed flux in an isolated copy, using the active OF mode.

Run after sourcing the ordinary or ADR OpenFOAM environment. This links the
remaining symbols from the installed HisaSolver; it is not a full solver build.
No installed library or OpenFOAM configuration is modified.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('output',type=Path,help='new directory for build and logs')
args=parser.parse_args()
root=Path(__file__).resolve().parents[2]
out=args.output.resolve()
out.mkdir(parents=True,exist_ok=False)
shutil.copytree(root/'src',out/'src',ignore=shutil.ignore_patterns('linux64*','lnInclude'))
for name in ('buildDefs.include','versionInfo.H'):
    shutil.copy2(root/name,out/name)
for path in out.rglob('*'):
    if path.is_file():
        path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
(out/'lib').mkdir()
(out/'src/Make/files').write_text('finiteVolume/fluxSchemes/ausmPlusUp/ausmPlusUpFluxScheme.C\n\nLIB = '+str(out/'lib/libHisaPressureFix')+'$(WM_AD_MODE)\n')
options=out/'src/Make/options'
options.write_text(options.read_text()+'\nLIB_LIBS += -L$(DAFOAM_ROOT_PATH)/OpenFOAM/sharedLibs -lHisaSolver$(WM_AD_MODE) -Wl,-z,defs\n')
with (out/'build.log').open('w') as log:
    result=subprocess.run(['wmake','-j2','libso',str(out/'src')],stdout=log,stderr=subprocess.STDOUT)
mode=os.environ.get('WM_AD_MODE','')
report={'mode':mode or 'native','returncode':result.returncode,'library':str(out/'lib'/f'libHisaPressureFix{mode}.so'), 'scope':'changed translation unit linked against installed HisaSolver; not a full CFD validation'}
(out/'build.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(result.returncode)
