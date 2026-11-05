#!/usr/bin/env python3
"""Serial pure-aero validation with explicit library provenance and raw gates."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--case',type=Path,required=True)
parser.add_argument('--library-dir',type=Path,required=True)
parser.add_argument('--alpha',type=float,default=2.)
parser.add_argument('--restart',action='store_true')
parser.add_argument('--no-wall-functions',action='store_true')
parser.add_argument('--adjoint',action='store_true')
parser.add_argument('--allow-unconverged-adjoint',action='store_true',help='Diagnostic only; does not mark the primal or gradient accepted')
parser.add_argument('--adjoint-restart',type=int,default=100)
parser.add_argument('--adjoint-full-lu',action='store_true',help='Small-case diagnostic: full residual FD Jacobian with LU preconditioning')
parser.add_argument('--probe',action='store_true',help='reload final_state.npz without primal iterations')
args=parser.parse_args()
case=args.case.resolve(strict=True)
libdir=args.library_dir.resolve(strict=True)
os.chdir(case)
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from dafoam import PYDAFOAM
from threadpoolctl import threadpool_limits
limits=threadpool_limits(limits=1,user_api='blas')
if MPI.COMM_WORLD.size!=1:
    raise RuntimeError('This validation script requires one MPI rank')
ap=AeroProblem(name='wing',mach=.85,altitude=10000.,alpha=args.alpha,
               areaRef=45.5,chordRef=3.25,evalFuncs=['cl','cd'])
velocity=[float(ap.V*np.cos(np.deg2rad(args.alpha))),0.,float(ap.V*np.sin(np.deg2rad(args.alpha)))]
options={'solverName':'DAHisaFoam','designSurfaces':['wing'],
 'primalMinResTol':1e-8,'primalMinResTolDiff':1.,'printInterval':20,
 'normalizeStates':{'U':float(ap.V),'p':float(ap.P),'T':float(ap.T),'nuTilda':float(ap.nu),'phi':1.},
 'primalBC':{'useWallFunction':not args.no_wall_functions},
 'checkMeshThreshold':{'maxAspectRatio':100000.,'maxNonOrth':90.,'maxSkewness':4.,'maxIncorrectlyOrientedFaces':0},
 'inputInfo':{'patchV':{'type':'patchVelocity','patches':['inout'],'flowAxis':'x','normalAxis':'z','components':['solver','function']}},
 'function':{},
 'adjEqnOption':{'gmresRestart':args.adjoint_restart,'gmresMaxIters':2000,'gmresRelTol':1e-8,'gmresAbsTol':1e-13,
                 'pcFillLevel':1,'jacMatReOrdering':'rcm','asmOverlap':1,
                 'hisaPCFlux':1,'hisaForceJSTFlux':0,'hisaRemoveTauMCOff':1},
 'adjPartDerivFDStep':{'State':1e-7}}
for field,value in [('U',velocity),('p',[float(ap.P)]),('T',[float(ap.T)]),('nuTilda',[float(ap.nu)])]:
    options['primalBC'][field+'0']={'variable':field,'patches':['inout'],'value':value}
for name,mode in [('cd','parallelToFlow'),('cl','normalToFlow')]:
    options['function'][name]={'type':'force','source':'patchToFace','patches':['wing'],
      'directionMode':mode,'patchVelocityInputName':'patchV','scale':float(1/(ap.q*ap.areaRef))}
if not args.restart:
    options['primalInitCondition']={'U':velocity,'p':float(ap.P),'T':float(ap.T),'nuTilda':float(ap.nu)}
(case/'options.json').write_text(json.dumps(options,indent=2)+'\n')
solver=PYDAFOAM(options=options,comm=MPI.COMM_WORLD)
loaded={}
for line in Path('/proc/self/maps').read_text().splitlines():
    if '/libHisaSolver' in line:
        path=Path(line.split()[-1]).resolve()
        if path.parent!=libdir:
            raise RuntimeError(f'Wrong HISA library loaded: {path}')
        loaded[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
if not any(p.endswith('libHisaSolver.so') for p in loaded) or not any(p.endswith('libHisaSolverADR.so') for p in loaded):
    raise RuntimeError(f'Both native and ADR patched libraries must be loaded: {loaded}')
solver.set_solver_input({'patchV':np.array([float(ap.V),args.alpha])})
initial_residual=None if args.probe else solver.getResiduals().copy()
if not args.probe:
    np.savez(case/'initial_state.npz',state=solver.getStates().copy())
start=time.monotonic()
if args.probe:
    solver.setStates(np.load(case/'final_state.npz')['state'])
else:
    solver()
elapsed=time.monotonic()-start
state=solver.getStates().copy(); residual=solver.getResiduals(); cells=len(state)//6
report={'alpha':args.alpha,'primal_fail':int(solver.primalFail),'raw_residual_max':float(np.max(np.abs(residual))),
 'raw_residual_l2':float(np.linalg.norm(residual)),'primal_seconds':elapsed,'cells':cells,'primal_iterations':0 if args.probe else int(round(solver.solver.getLatestTime())),
 'temperature_min':float(state[4*cells:5*cells].min()),'temperature_max':float(state[4*cells:5*cells].max()),
 'functions':{name:float(solver.solver.calcFunction(name)) for name in ('cd','cl')},'loaded_libraries':loaded,
 'loaded_dafoam_libraries':{str(Path(line.split()[-1]).resolve()):hashlib.sha256(Path(line.split()[-1]).read_bytes()).hexdigest() for line in Path('/proc/self/maps').read_text().splitlines() if '/libDASolver' in line},
 'accepted':bool(solver.primalFail==0 and np.all(np.isfinite(residual)) and np.max(np.abs(residual))<1e-8),
 'gradient_verified':False,'initialization':'saved fields' if args.restart else 'uniform freestream',
 'adjoint_settings_configured_but_not_run':options['adjEqnOption']}
if initial_residual is not None:
    np.savez(case/'initial_residual.npz',residual=initial_residual)
    report['initial_raw_residual_l2']=float(np.linalg.norm(initial_residual))
    report['fixed_definition_relative_l2']=float(np.linalg.norm(residual)/max(np.linalg.norm(initial_residual),1e-300))
np.savez(case/'final_state.npz',state=state,residual=residual)
(case/('probe_report.json' if args.probe else 'primal_report.json')).write_text(json.dumps(report,indent=2)+'\n')
print('VALIDATION_REPORT '+json.dumps(report),flush=True)

if args.adjoint and not report['accepted'] and not args.allow_unconverged_adjoint:
    (case/'adjoint_report.json').write_text(json.dumps({'executed':False,'reason':'Base primal failed raw residual gate','gradient_verified':False},indent=2)+'\n')
if args.adjoint and (report['accepted'] or args.allow_unconverged_adjoint):
    from petsc4py import PETSc
    solver.setStates(state)
    patch_velocity=np.array([float(ap.V),args.alpha])
    solver.set_solver_input({'patchV':patch_velocity})
    solver.solver.runColoring()
    pc_matrix=PETSc.Mat().create(MPI.COMM_WORLD)
    solver.solver.calcdRdWT(0 if args.adjoint_full_lu else 1,pc_matrix)
    solver.solverAD.initializedRdWTMatrixFree()
    ksp=PETSc.KSP().create(MPI.COMM_WORLD)
    solver.solverAD.createMLRKSPMatrixFree(pc_matrix,ksp)
    if args.adjoint_full_lu:
        ksp.getPC().setType('lu')
        if PETSc.Sys.hasExternalPackage('mumps'):
            ksp.getPC().setFactorSolverType('mumps')
    adj={'executed':True,'unconverged_base_diagnostic':not report['accepted'],'gmres_restart':options['adjEqnOption']['gmresRestart'],
         'pc_matrix':'calcdRdWT(isPC=0): full-residual colored FD Jacobian transpose' if args.adjoint_full_lu else 'calcdRdWT(isPC=1): colored finite-difference approximate transposed residual Jacobian',
         'pc_shape':pc_matrix.getSize(),'pc_nonzeros':pc_matrix.getInfo()['nz_used'],
         'preconditioning':'right LU of full-residual FD matrix (small-case diagnostic)' if args.adjoint_full_lu else 'right ASM(overlap=1), local ILU(1), RCM ordering',
         'pc_residual':'same high-order residual as AD operator' if args.adjoint_full_lu else 'first-order LaxFriedrichs convection; off-diagonal tauMC removed; hisaForceJSTFlux=0',
         'gradient_verified':False,'functions':{}}
    for name in ('cd','cl'):
        rhs=np.zeros(len(state))
        solver.solverAD.calcJacTVecProduct('states','stateVar',state,name,'function',np.ones(1),rhs)
        rhs_vec=solver.array2Vec(rhs); psi=rhs_vec.duplicate(); psi.set(0.)
        solve_start=time.monotonic()
        fail=solver.solverAD.solveLinearEqn(ksp,rhs_vec,psi)
        solve_seconds=time.monotonic()-solve_start
        seed=psi.getArray().copy(); product=np.zeros(len(state))
        solver.solverAD.calcJacTVecProduct('states','stateVar',state,'residuals','residual',seed,product)
        explicit=np.zeros(2);implicit=np.zeros(2)
        solver.solverAD.calcJacTVecProduct('patchV','patchVelocity',patch_velocity,name,'function',np.ones(1),explicit)
        solver.solverAD.calcJacTVecProduct('patchV','patchVelocity',patch_velocity,'residuals','residual',seed,implicit)
        adj['functions'][name]={'d_dalpha':float((explicit-implicit)[1]),'linear_fail':int(fail),'solve_seconds':solve_seconds,
          'gmres_iterations':ksp.getIterationNumber(),'ksp_reason':int(ksp.getConvergedReason()),
          'true_linear_relative_residual':float(np.linalg.norm(product-rhs)/max(np.linalg.norm(rhs),1e-30))}
        (case/'adjoint_report.json').write_text(json.dumps(adj,indent=2)+'\n')
        np.savez(case/('adjoint_'+name+'.npz'),psi=seed,rhs=rhs,product=product)
        rhs_vec.destroy();psi.destroy()
    ksp.destroy();pc_matrix.destroy()
