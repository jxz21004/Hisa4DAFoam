# L3/L2 pure-aerodynamic convergence and AoA-gradient validation

This is a serial, isolated validation for the user's wing meshes, not an
optimization or aerostructural run. It requires the existing workspace's
`dafoam_case/converter_package_validation/{L3,L2}_final` meshes and level-specific
field templates documented in `prepare.py`. CGNS inputs and the installed solver
are not modified. Reports record the paths and SHA256 of the libraries actually
mapped into each Python process.

Prepare a NEW work directory, then build both complete HISA libraries:

```bash
python tests/wingAero/prepare.py /absolute/path/to/new-validation
# Ordinary OpenFOAM environment
wmake -j2 libso /absolute/path/to/new-validation/build/src
# Separate shell, reverse-AD OpenFOAM environment
source "$DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-AD/etc/bashrc" WM_AD_MODE=ADR
wmake -j2 libso /absolute/path/to/new-validation/build/src
```

Return to the ordinary DAFoam environment, then run:

```bash
python tests/wingAero/run_suite.py /absolute/path/to/new-validation
```

`run_suite.py` prepends the isolated library directory to each child process's
LD_LIBRARY_PATH. It runs L3 and L2 in separate serial processes, starts from
uniform freestream, and requires maximum absolute raw residual < 1e-8 before
running converged ±0.001 degree AoA cases and adjoints. Function directions also
vary with AoA. A successful native exit is not treated as proof of convergence.

The primal GMRES uses nKrylov=8, maxIter=60 restart cycles, LUSGS, LaxFriedrichs
inviscid and laplacian viscous Jacobians. SA uses its separate smoothSolver.
Adjoints, when the prerequisite flows converge, use the actual reverse-AD Jᵀ,
PETSc GMRES restart=100/maxIters=2000, and the colored finite-difference
`calcdRdWT(isPC=1)` matrix: first-order LaxFriedrichs PC flux, off-diagonal tauMC
removed, right ASM(1)/ILU(1), RCM ordering. `hisaForceJSTFlux=0` is explicit.
The script saves actual linear iteration counts and independently evaluates the
true linear residual. Configured-but-unused settings are not executed results.

The completed run also contains a separate L3 restart comparison: clear
`0/pseudoCoField[.gz]` and `0/uniform`, set pseudoCoNumMax=5/nKrylov=32/maxIter=10,
and run another 2000 steps with `run_case.py --restart`. Its directory is
`L3_lowCFL_clean/base`. `summarize.py WORK` summarizes all three completed runs
and plots the maximum component of HiSA's initial normalized linear residual;
that plotted residual is distinct from the raw DAFoam residual acceptance gate.

Completed ac05bed results: [REPORT.md](results/REPORT.md),
[convergence plot](results/convergence.png), and JSON/configuration snapshots in
`results/`. All three selected runs failed the steady residual gate; no steady
AoA gradient accuracy result is claimed. Full logs and fields remain in
`dafoam_case/validation/hisa_ac05bed` outside the repository.
Configuration snapshots normalize line endings and trailing whitespace; mesh
hashes refer to the unchanged actual input files.
