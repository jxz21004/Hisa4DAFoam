#!/usr/bin/env python3
"""Instrument an isolated existing build; production sources stay untouched."""
import argparse,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('base',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
shutil.copytree(a.base,a.output,symlinks=True)
f=a.output/'src/Make/files';f.write_text(f.read_text().replace(str(a.base.resolve()),str(a.output.resolve())))
f=a.output/'src/finiteVolume/fluxSchemes/ausmPlusUp/ausmPlusUpFluxScheme.C';s=f.read_text()
s=s.replace('#include "cellFaceFunctions.H"','#include "cellFaceFunctions.H"\n#include "linear.H"\n#include <fstream>\n#include <cstdlib>')
needle='    Up = U_f()*mesh_.magSf();'
assert s.count(needle)==1
s=s.replace('surfaceFieldSelect(U_L, U_R, Mach_1_2(), 0)', 'surfaceFieldSelect(U_L(), U_R(), Mach_1_2(), 0)')
s=s.replace(needle,needle+'''
#ifndef CODI_ADR
    // Diagnostic-only experiment, never installed or enabled in production.
    if (std::getenv("HISA_PROBE_LINEAR_VISCOUS"))
    {
        Up = linearInterpolate(U_)*mesh_.magSf();
    }
    if (const char* path = std::getenv("HISA_PROBE_CSV"))
    {
        std::ofstream out(path);
        out.precision(17);
        out << "face,owner,neighbour,mach,area,massFlux,energyFlux,pressure,ULx,ULy,ULz,URx,URy,URz,Upx,Upy,Upz\\n";
        forAll(mesh_.neighbour(), face)
        {
            out << face << ',' << mesh_.owner()[face] << ',' << mesh_.neighbour()[face]
                << ',' << Mach_1_2()[face] << ',' << mesh_.magSf()[face]
                << ',' << rhoa_LR[face]*mesh_.magSf()[face]
                << ',' << rhoah_LR[face]*mesh_.magSf()[face] << ',' << p_1_2[face];
            for (label d=0; d<3; ++d) out << ',' << U_L()[face][d];
            for (label d=0; d<3; ++d) out << ',' << U_R()[face][d];
            for (label d=0; d<3; ++d) out << ',' << Up[face][d];
            out << '\\n';
        }
    }
#endif
''')
f.write_text(s)
