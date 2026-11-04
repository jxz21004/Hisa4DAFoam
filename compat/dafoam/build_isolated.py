#!/usr/bin/env python3
"""Copy this installation's cached build, apply the SA patch, build outside it."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, default=Path('/home/dafoamuser/dafoam'))
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
root, out = a.root.resolve(), a.output.resolve()
# Require a fresh destination; never merge with or overwrite an installation.
out.mkdir(parents=True, exist_ok=False)
for name in ('adjoint', 'include'):
    shutil.copytree(root / 'repos/dafoam/src' / name, out / name, symlinks=True)
patch = Path(__file__).with_name('sa_hisa_second_ring_velocity.patch')
# Patch paths start with src/; the copied tree starts at adjoint/.
subprocess.run(['patch', '-p2', '--forward', '-i', str(patch.resolve())], cwd=out, check=True)
(out / 'lib').mkdir()
files = out / 'adjoint/Make/files'
lines = files.read_text().splitlines()
assert sum(line.strip().startswith('LIB =') for line in lines) == 1
files.write_text('\n'.join('LIB = ' + str(out / 'lib/libDASolver$(WM_AD_MODE)')
                          if line.strip().startswith('LIB =') else line for line in lines) + '\n')
env = dict(os.environ, DAFOAM_ROOT_PATH=str(root))
for mode, tree in (('', 'OpenFOAM-v2506'), ('ADR', 'OpenFOAM-AD')):
    bashrc = root / 'OpenFOAM' / tree / 'etc/bashrc'
    with (out / ('build_' + (mode or 'native') + '.log')).open('w') as log:
        subprocess.run(['bash', '-c', 'source "$1" "WM_AD_MODE=$2"; wmake -j2 libso "$3"',
                        'build', str(bashrc), mode, str(out / 'adjoint')],
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
