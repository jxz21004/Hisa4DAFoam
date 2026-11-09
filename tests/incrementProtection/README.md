# CFL3D-style conservative increment protection

Opt-in pseudoTime settings: `incrementProtection true; incrementProtectionReportInterval 100;`.

Uses CFL3D update.F's scalar coefficient: 1 for relative increment r >= -0.2, otherwise 1/(0.6-2*r). Normal updates are unchanged; negative relative changes approach -0.5. No exponentials/logarithms. This is continuous but not C1.

Adaptation: apply min(alpha_density,alpha_linearized_pressure) to ALL conservative increments. Check actual density and internal-energy density remain >=50% of the old state; halve alpha only if necessary. Reject nonfinite proposals. This is not a literal primitive-variable update from CFL3D, and does not guarantee convergence or limit cumulative heating/acceleration. Roe-only for now (its flux validates the perfect-gas energy assumptions).

Old separate local density/internal-energy clipping is bypassed in this mode. Thermodynamic emergency bounds remain. No residual/Jacobian/flux re-evaluation, no per-step MPI reductions introduced: diagnostics/reductions only at the configured interval (0 disables). Uses one cell loop and local arithmetic. Disabled by default. Do not combine with matrixFree or limitedUpdate.

Run kernel checks:
```
g++ -std=c++17 -O2 -I src/modules/hisa tests/incrementProtection/test_kernel.cpp -o /tmp/test_increment
/tmp/test_increment
```
Native validation only; AD compilation/gradient validation not yet performed. The protection belongs to the primal update, not the DAFoam steady residual tape.

Reference: https://github.com/nasa/CFL3D/blob/master/source/cfl3d/libs/update.F
