# discreteLaxFriedrichs orientation compatibility

OpenFOAM-v2506: geometric mesh weights and reconstructed interpolation weights
can inherit an oriented flag. Treating these scalar interpolation coefficients
as oriented causes the implicit gradient/divergence operators to add fields
with incompatible orientation flags and abort before the first update.

The fix sets w, wu, we, wrho unoriented, guarded by OPENFOAM >= 1712, matching
the convention already used for w in laxFriedrichs. No coefficient values,
explicit residual scheme, turbulence model, or other Jacobian branch changes.

Validation:
- Native and ADR HiSA libraries rebuilt in an isolated workspace directory.
- Before: discrete+CFL5 aborts on the first iteration with an oriented/unoriented
  addition error. Clearing only w gets past the gradient but fails in divMult;
  the reconstructed weights must be corrected too.
- After: discrete+CFL5 with separate central LUSGS matrix completes 300 steps
  from the same L3 state, residual L2 80068.8245 -> 31539.5806.
- Native/ADR residual parity relative error: 1.2097e-16.
- CFL100 remains numerically unstable. This fix is not a nonlinear robustness
  or convergence fix; it enables testing the previously unusable branch.

Reproduction scripts and numerical results are in the sibling DAFoam repo:
`tests/hisaConnectivity/operatorSweep/README.md`, `operator_sweep.py`, and
`cfl_sweep.py`. Use --jacobian discreteLaxFriedrichs and
--pc-jacobian centralLaxFriedrichs. Libraries used for after-fix tests:
`dafoam_case/validation/hisa_discrete_fix/build/lib` under the workspace root.
The case reports retain native/ADR hashes; system installation is untouched.
