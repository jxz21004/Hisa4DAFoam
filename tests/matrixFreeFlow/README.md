# Matrix-free native flow solve (experimental, default off)

`flowSolver/GMRES/matrixFree true` enables residual-based products while keeping
the assembled matrix for LUSGS and the legacy outer residual normalization.
See sibling DAFoam `tests/hisaConnectivity/matrixFreeFlow/README.md` and
`results.json` for exact configuration, scope, validation and unsuccessful
convergence results. This is a frozen-turbulence flow experiment, not a fully
coupled NK solver or a production default.

Implementation:
- hisaSolver provides an optional matrix-product callback; no callback uses
  the original matrix multiplication.
- GMRES uses the callback for Arnoldi and actual linear-residual checks.
- matrixFreeSolve.H samples the native high-order conservative residual with
  central finite differences; its step is scaled by global state/direction
  magnitudes. Pseudo-time coefficients and turbulence are frozen.
- Trial residual evaluation reuses native field/flux/residual updates and
  restores the baseline. A zero-step reproducibility gate aborts on drift.
- Optional half-step comparison checks a first Krylov direction per outer step.
- The candidate flow increment is backtracked against a scaled steady residual
  norm; failure restores the initial flow state. SA remains a separate update.

Native/ADR builds passed. Turning the option off reproduces old-library state
and residual arrays exactly over 20 iterations. Three directional checks had
maximum relative step-halving difference 3.72e-8 and zero-step drift 9.63e-13.
At large CFL, the true linear problems did not reach tolerance within 120
Krylov steps. All flow tests remain unconverged; do not interpret primalFail=0
as convergence. The stored reports include rejected steps and evaluation cost.

Tests used the existing copied build at
`dafoam_case/validation/hisa_matrix_free/build` under the workspace. The changed
headers require rebuilding the HiSA solver templates and hisaModule for both
native and ADR; `src/Make/files` in that build targets its isolated lib directory.
DAFoam's existing Python wrappers load those libraries via LD_LIBRARY_PATH.
