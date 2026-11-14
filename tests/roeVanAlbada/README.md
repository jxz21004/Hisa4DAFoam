# Smooth squared-slope van Albada type limiter

Optional `roeLimiter vanAlbadaType` for Roe UMUSCL. This follows the squared
slope modification of the regularized van Albada factor described by
Engels-Putzka, Backhaus and Frey (2014), equation13 applied to squared slopes:
https://elib.dlr.de/90138/1/Engels-Putzka_2014_ECCOMAS.pdf
It is NOT a claim to reproduce elsA's implementation or the sign-clipped TVD
van Albada limiter.

For owner cell i and neighbor j, define D=Cj-Ci, b=qj-qi, a=2*grad(q)i.D-b.
These are directional backward/forward differences on a virtual stencil.
The existing geometrically corrected UMUSCL face increment delta is multiplied by

    phi = (2*a^2*b^2 + e^2) / (a^4 + b^4 + e^2)
    e = q_ref^2 * K^3 * V / L_ref^3 > 0
    q_face = qi + phi*delta

The opposite reconstruction uses its own gradient, volume and reversed D.
No division by slope, sign switch, neighbor extrema or hard min/max clipping is
used in this limiter. The rational factor is smooth for e>0 and bounded by1.
Equal-magnitude opposite slopes also give phi=1: smooth extrema are retained,
but strict TVD and neighbor bound preservation are deliberately not guaranteed.
Roe positivity fallback remains independent and unchanged; other solver branches
(including entropy corrections) can still be nonsmooth. The approximate primal
Jacobian remains the existing frozen Roe approximation.

Example with dimensional RAE2822 reference scales:

```
roeLimiter vanAlbadaType;
roeVanAlbadaCoeff 1;
roeVanAlbadaReferenceLength 1;
roeVanAlbadaReferenceRho 0.4309079200663052;
roeVanAlbadaReferenceU 220;
roeVanAlbadaReferenceP 27615;
```

The five numeric values default to1 and must be positive. The coefficient is K,
not a CFL factor. Larger K weakens limiting. Fixed reference scales are shared
by base/perturbed cases; they do not track solution extrema. The local volume
is active for shape derivatives and exchanged at processor boundaries. Because
V is a 3D cell volume, quasi-2D extrusion thickness affects regularization.
Legacy default remains hard limiting. Existing roeVenkat* settings affect only
Venkatakrishnan; this option has independent roeVanAlbada* settings.

Kernel tests (CoDiPack include path required):

```
g++ -std=c++17 -O2 -I src/finiteVolume/fluxSchemes/roe -I "$CODI_INCLUDE" tests/roeVanAlbada/test_kernel.cpp -o /tmp/test_roe_van_albada
/tmp/test_roe_van_albada
```

Covers coefficient bounds, symmetric/zero/opposite slopes, dimensional scaling,
linear-exact reconstruction including non-midpoint faces, and CoDi derivatives
of delta, both slopes and the geometric regularizer. Project parallel and AD
studies live under validation/hisa_adflow_comparison/roe_van_albada.

## Ordinary signed-slope option

`roeLimiter vanAlbada` uses the same virtual directional stencil and existing
UMUSCL increment, but the factor is `(2*a*b+e)/(a*a+b*b+e)` (no squared-slope
substitution). This is the regularized symmetric van Albada factor in equation13
of the cited reference, not the additional squared-slope modification. The
existing roeVanAlbada* scales apply. No max(0,phi) is applied; phi can be negative
for opposite slopes (smooth but not strict TVD). It is not claimed to exactly
reproduce elsA. The stopped squared-type flow run must not be relabeled as this
option. test_ordinary.cpp verifies signed factors, ratio2->0.8 and3->0.6,
linear reconstruction, scaling, zero slopes and all CoDi input derivatives.
Use the same compile command as above with test_ordinary.cpp.
