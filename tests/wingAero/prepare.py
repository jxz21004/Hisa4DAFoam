#!/usr/bin/env python3
"""Prepare an isolated full-library build and L3/L2 cases in this workspace."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('work',type=Path)
args=parser.parse_args()
repo=Path(__file__).resolve().parents[2];project=repo.parent
work=args.work.resolve();work.mkdir(parents=True,exist_ok=False)
build=work/'build';build.mkdir()
shutil.copytree(repo/'src',build/'src',ignore=shutil.ignore_patterns('linux64*','lnInclude'))
for name in ('buildDefs.include','versionInfo.H'):shutil.copy2(repo/name,build/name)
for path in build.rglob('*'):
    if path.is_file():path.write_bytes(path.read_bytes().replace(b'\r\n',b'\n'))
(build/'lib').mkdir()
p=build/'src/Make/files';p.write_text(re.sub(r'LIB = .*','LIB = '+str(build/'lib/libHisaSolver')+'$(WM_AD_MODE)',p.read_text()))
meshes={}
for level in ('L3','L2'):
    template=project/'dafoam_case/validation'/('L3_convert_margin01' if level=='L3' else 'L2_basic_block')
    case=work/level/'base';case.mkdir(parents=True)
    for folder in ('system','constant','0'):shutil.copytree(template/folder,case/folder)
    shutil.rmtree(case/'constant/polyMesh')
    shutil.copytree(project/'dafoam_case/converter_package_validation'/(level+'_final')/'constant/polyMesh',case/'constant/polyMesh')
    if (case/'0/uniform').exists():shutil.rmtree(case/'0/uniform')
    for name in ('pseudoCoField','pseudoCoField.gz'):
        if (case/'0'/name).exists():(case/'0'/name).unlink()
    p=case/'system/fvSchemes';s=p.read_text().replace('AUSMPlusUpContinuous','AUSMPlusUp')
    s=re.sub(r'(reconstruct\((?:rho|U|T)\))\s+[^;]+;',r'\1 wVanLeer;',s)
    s=re.sub(r'viscousFaceVelocity\s+linear;','',s);p.write_text(s)
    p=case/'system/controlDict';s=p.read_text()
    for key,value in [('endTime',2000),('writeInterval',2000)]:s=re.sub(rf'\b{key}\s+[^;]+;',f'{key} {value};',s)
    s=re.sub(r'\bHisaPrintInterval\s+[^;]+;','',s);s+='\nHisaPrintInterval 1;\n';p.write_text(s)
    p=case/'system/fvSolution';s=p.read_text()
    for key,value in [('pseudoCoNum',1),('pseudoCoNumMin',.1),('pseudoCoNumMax',100),('nKrylov',8),('maxIter',60)]:
        s=re.sub(rf'\b{key}\s+[^;]+;',f'{key} {value};',s)
    p.write_text(s)
    meshes[level]={n:hashlib.sha256((case/'constant/polyMesh'/n).read_bytes()).hexdigest() for n in ('points','faces','owner','neighbour','boundary')}
(work/'mesh_hashes.json').write_text(json.dumps(meshes,indent=2)+'\n')
print('Prepared',work)
