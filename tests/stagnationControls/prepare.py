#!/usr/bin/env python3
"""Prepare single-factor L3 controls and a diagnostic zero-slope reconstruction."""
from pathlib import Path
import shutil,re,json
root=Path(__file__).resolve().parents[3];w=root/'dafoam_case/validation/hisa_stagnation_controls';w.mkdir()
base=root/'dafoam_case/validation/hisa_viscous_work_fix/build';build=w/'build';shutil.copytree(base,build,symlinks=True)
f=build/'src/Make/files';f.write_text(f.read_text().replace(str(base),str(build)))
f=build/'src/finiteVolume/reconstructionSchemes/reconstructionScheme/reconstructionSchemes.C';s=f.read_text();needle='makeLimitedReconstructionScheme(Minmod, MinmodLimiter)'
s=s.replace(needle,needle+'''

// Diagnostic first-order reconstruction: zero limiter gives owner/neighbour
// face states through the existing reconstruction and boundary machinery.
template<class LimiterFunc>
class controlFirstOrderLimiter : public LimiterFunc
{
public:
    controlFirstOrderLimiter(Istream&) {}
    scalar limiter
    (
        const scalar, const scalar,
        const typename LimiterFunc::phiType&, const typename LimiterFunc::phiType&,
        const typename LimiterFunc::gradPhiType&, const typename LimiterFunc::gradPhiType&,
        const vector&
    ) const { return scalar(0); }
};
makeLimitedReconstructionScheme(controlFirstOrder, controlFirstOrderLimiter)
''');f.write_text(s)
original=root/'dafoam_case/validation/hisa_viscous_fix_L3_2000/L3/base'
for name in ['first_order','high_order_cfl5']:
 c=w/name/'L3/base';c.mkdir(parents=True);(w/name/'build').symlink_to(build if name=='first_order' else base,target_is_directory=True)
 for d in ['0','constant','system']:shutil.copytree(original/d,c/d)
 if name=='first_order':
  f=c/'system/fvSchemes';f.write_text(re.sub(r'(reconstruct\((?:rho|U|T)\))\s+wVanLeer;',r'\1 controlFirstOrder;',f.read_text()))
 else:
  f=c/'system/fvSolution';f.write_text(re.sub(r'\bpseudoCoNumMax\s+[^;]+;','pseudoCoNumMax 5;',f.read_text()))
(w/'settings.json').write_text(json.dumps(dict(base_commit='6c1d633',iterations=2000,initialization='uniform freestream',rtol=.05,atol=1e-12,nKrylov=40,controls=dict(first_order='Only reconstruct(rho/U/T) replaced by zero-slope reconstruction',high_order_cfl5='Only pseudoCoNumMax changed from 100000 to 5')),indent=2)+'\n')
print(w)
