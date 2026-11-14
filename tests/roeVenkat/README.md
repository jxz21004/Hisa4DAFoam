# Roe UMUSCL Venkatakrishnan reconstruction

Select in `system/fvSchemes`:

```
fluxScheme Roe;
roeLimiter venkatakrishnan; // none | hard | venkatakrishnan
roeVenkatCoeff 1;
roeVenkatReferenceLength 1; // metres
roeVenkatReferenceRho 1;   // kg/m^3
roeVenkatReferenceU 1;     // m/s, common to all velocity components
roeVenkatReferenceP 1;     // Pa
```

All five numeric options default to 1 and must be positive when enabled.
Use fixed characteristic scales for a dimensional case (usually chord and
freestream rho, speed and pressure); do not recompute scales from current
solution extrema. This implementation is facewise Venkatakrishnan applied to
existing UMUSCL increments, not a cell-wide minimum of face limiter factors.

Let d be the unlimited face-minus-cell increment and a the upper/lower neighbor
bound-minus-cell value selected by the sign of d. The limited increment is

    d * (a*a + 2*a*d + eps2) / (a*a + a*d + 2*d*d + eps2)
    eps2 = q_ref^2 * K^3 * V / L_ref^3

This is the classical rational limiter, written without division by d.
It can allow small overshoots and coefficients above one; no extra hard clipping
is applied. Larger K relaxes limiting. V is the local 3D cell volume, including
extrusion thickness for quasi-2D meshes: changing thickness changes smoothing.
The formula follows the Venkatakrishnan formulation discussed in
https://web.stanford.edu/group/ctr/ResBriefs09/03_pecnik.pdf ; dimensional state
and length scales are explicit here. This does not guarantee global C1 residuals:
neighbor extrema selection and other solver branches remain.

Missing roeLimiter preserves legacy roeLimit true -> hard, false -> none.
If both keys are explicitly supplied they must agree on whether limiting is on.
Internal and processor faces share the same formula. The neighbor smoothing
scale is exchanged with processor patches; geometry remains active for AD.
Physical boundary fluxes, Roe entropy treatment, approximate primal Jacobian,
and positivity fallback are unchanged. The native and ADR library paths both
use this implementation. A converged old-limiter field is not a converged
solution of the new residual and must be reconverged for gradient verification.

Kernel regression (set CODI_INCLUDE to the CoDiPack include directory):

```
g++ -std=c++17 -O2 -I src/finiteVolume/fluxSchemes/roe -I "$CODI_INCLUDE" tests/roeVenkat/test_kernel.cpp -o /tmp/test_roe_venkat
/tmp/test_roe_venkat
```

Tests cover known values, zero-increment behavior, sign and dimensional scaling,
and CoDi forward derivatives of all inputs, including the smoothing parameter.
Project frozen-field/parallel studies live outside this repository under
validation/hisa_adflow_comparison/roe_venkat.
