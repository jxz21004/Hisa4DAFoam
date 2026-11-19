# DAFoam HISA hybrid primal solver

`HisaPetscNK.H` owns the passive C++ PETSc NK algorithm. DAFoam's
`DAHisaFoamHybrid.C` supplies state access, the full flow/SA residual, and the
sparse approximate preconditioner. Python configures and launches the solve;
there is no Python callback in the NK iteration. This is an opt-in experimental
integration, not the standalone HISA executable's NK implementation.

The current adapter supports steady RANS with U, p, T and nuTilda, state ordering
`state`, and positive `normalizeStates` for all four variables. Use the native
primal instance; reverse AD remains available for residual/adjoint evaluation.
It does not differentiate through the nonlinear iteration.

Add to `system/fvSolution`:

```
hisaHybrid
{
    enabled true;
    nkSwitchTol 1e-6;
    ptcRecoverySteps 10;
    steadyRelTol 1e-8;
    steadyAbsTol 1e-14;
    ncycles 8000;
    ptcMinIterations 10;
    nkSubspace 60;
    nkMaxLinearIterations 60;
    nkASMOverlap 1;
    nkILUFill 2;
    nkJacobianLag 20;
    nkLinearTol 0.3;
    nkUseEW true;
    monitorInterval 10;
    monitors (relative absolute cfl ksp linear alpha);
    detailed false;
}
```

These are the defaults except `enabled` (false when absent). RCM is explicitly
selected inside each ASM subdomain. NK uses right-preconditioned GMRES, WP
finite-difference Jv, an EW tolerance capped at 0.8, and cubic Armijo backtracking.
The sparse DAFoam PC is transposed into primal orientation and row-scaled; it is
not the exact residual Jacobian. Preconditioners are rebuilt on stage entry and
restart, then reused for `nkJacobianLag` updates. Rejected line-search steps
restore the base state and return to PTC, rather than accepting a worsening
minimum step. This intentionally differs from the earlier experimental script.

`absolute` is the global L2 norm of the fixed-scaled full residual, identical
in definition to the NK line-search merit residual. Each entry is scaled by the
inverse initial RMS of its U/p/T/nuTilda block (floor 1e-30). `relative` divides
this total norm by its original total norm. Final convergence requires total
relative <= steadyRelTol OR total absolute <= steadyAbsTol. Individual blocks
are diagnostic only. `blocks` prints four unscaled DAFoam block RMS values;
`combined` is an alias for total relative. `absolute` is a scaled discrete norm,
not a dimensional conservation residual or the native GMRES inner residual.
Stage-specific monitors: `cfl` is printed only for PTC. In PTC, `ksp`
prints actual native GMRES Krylov steps and restart cycles; `linear` prints
`linRelMax` (maximum component-wise final/initial linear residual ratio) and
`linearStop` (`absolute`, `relative`, or `maxIter`). These reuse the solver's
existing residual checks. In NK, `ksp`, `linear`, and `alpha` print PETSc steps,
recomputed linear defect ratio, and accepted line-search step respectively.
PTC never prints `alpha`; NK never prints placeholder CFL. Solve statistics
are omitted before the first solve and when a stage change makes them stale.
PTC `linRelMax` and NK `linRel` use different norms and should not be compared
numerically as the same metric.

`ncycles` limits cumulative PTC + NK outer attempts, including across restarts.
There is no separate NK outer limit. `ptcMinIterations` (default10) requires at
least this many cumulative PTC updates before NK entry, including when restarting
an NK checkpoint with a newly increased PTC minimum. It does not prevent stopping
if the final convergence tolerance is already met. `ptcRecoverySteps` applies
separately after every NK fallback. `nkMaxLinearIterations` still limits each
inner GMRES solve. Old maxOuterIterations/nkMaxIterations keys are no longer used.

PTC uses native HISA flow iteration and DAFoam SA correction. In hybrid mode,
its SER target is now initialCFL / totalRelativeResidual, using exactly the fixed
scaled full residual and baseline used for NK switching and stopping. Existing
pseudoCoNumMin/Max and per-step increase/decrease factors still apply; local
CFL damping/increment protection remain unchanged. No extra CFL reduction is
applied solely because NK failed. Initial CFL is checkpointed with the baseline;
legacy checkpoints use configured pseudoCoNum on their first migration. The
current pseudoCoNum/pseudoCoField continue to use native checkpoint storage.
Use resetPseudo=false and limitedUpdate=false: competing CFL controllers are
rejected explicitly. Native non-hybrid HISA retains its original SER rule. NK does not use that CFL as a step
size: it uses line-search alpha. If the total relative residual exceeds `nkSwitchTol`, or NK
rejects its step, the next update uses PTC. Recovery steps prevent immediate
chatter. The same `nkSwitchTol` controls entry (ratio <= threshold) and return
(ratio > threshold); there is no separate return tolerance.

Each field checkpoint includes `<time>/uniform/hisaHybridState` with the fixed
block scales, initial total residual, stage and cumulative outer/NK counts. Preserve it on every MPI
rank along with fields. Restart with the same mesh and partition; use
`startTime`/`latestTime` normally and do not force decomposition over saved
processor fields. The adapter refuses a nonzero-time restart without a baseline
unless `resetBaseline true` explicitly requests a new solve. It checks global
cell count and scales, but does not fingerprint mesh topology. Changing the
physical problem requires a new baseline. PC factors and EW history are rebuilt,
so restart is not expected to reproduce the uninterrupted trajectory bitwise.

Iteration limits are cumulative, including rejected NK attempts. `endTime` is
also respected. Reaching a limit returns primal failure/status rather than
claiming convergence. Default console output is every ten outer iterations;
stage changes and final status always print. `detailed true` restores native
solver chatter. OpenFOAM function-object file output intervals remain separately
configured in controlDict. Disable duplicate function-object console monitors
if unwanted. With `enabled false`, the original solvePrimal path is used.

Legacy hybrid checkpoints without initialTotalResidual reconstruct it from their
persisted original block RMS and global cell count, never from the restart field.
PTC counts are derived from cumulative outer minus NK attempts.

Monitor `cfl` (enabled by default) prints the PTC global pseudoCoNum and local
pseudoCoField min/max at the reported state. NK prints CFL n/a and uses alpha.

Linear stagnation guard: `nkLinearStagnationTol` defaults to0.99 (valid range
0 < value <= 1). If the independently recomputed linear defect ratio reaches
this value or is nonfinite, restore the physical base state, rebuild the PC,
and retry the same linear problem once with the same EW tolerance. If still
stagnant, return a rejected update to the adapter, which falls back to PTC.
No line search is run on that stagnant direction. Retry counts as one outer
attempt; `KSP`/`totalKSP` includes both linear solves (up to twice the configured
linear limit). An immediate event records first/retry linRel; normal successful
steps retain the existing PC lag. No extra CFL reduction is introduced.


Hybrid PTC SER tuning in `pseudoTime` (independent of initial CFL):

```
pseudoCoSafetyFactor      1;
pseudoCoResidualExponent  1;
pseudoCoLogDamping        1;
```

Target = safetyFactor * savedInitialCFL / totalRelative^exponent. Both default
1 preserve the previous rule; safety must be finite >0, exponent finite >=0.
Exponent0 gives a constant target. Existing per-step limits and CFL bounds
still apply, including at startup. The initial CFL field is still initialized
from pseudoCoNum; a low target can then reduce it subject to the lower bound.
A safety factor0.1 turns an otherwise50 target into5 when exponent1. Exponent
changes how strongly accumulated residual reduction grows the target. These
settings affect hybrid PTC only, not NK line search or non-hybrid native SER.
They are user controls read from the dictionary, not saved baseline quantities.
Changing them on restart intentionally changes subsequent targets. Default
exponent1 avoids pow; otherwise one scalar power per PTC iteration per rank
is used (not one per cell/face). No wall-time speedup/slowdown claim is made;
changes to the nonlinear iteration count may dominate their arithmetic cost.

Optional `pseudoCoLogDamping` omega is finite in (0,1], default1. Before
existing rate limits and final CFL bounds, replace target/current by
(target/current)^omega. This is log-CFL interpolation, not a new residual
evaluation or a PD controller. Smaller omega gives stronger damping.
Omega1 skips the extra pow and preserves old arithmetic. The current saved
CFL is sufficient for restart; no additional history is required. This applies
to hybrid PTC only. Example trial: omega0.2, safety0.5, exponent0.5.

Optional stagnation integral: `pseudoCoStagnationGain` defaults0 (disabled),
`pseudoCoStagnationDecay` defaults0.002 (positive log decrease per PTC call).
For consecutive PTC residual evaluations s=log(previousRelative/relative),
I=max(0,I+decay-s). Multiply the SER target by exp(-gain*I), then apply
log damping, rate limits and absolute bounds. At minimum CFL, positive
integral increments are frozen; negative increments can unwind. First call
without history only seeds previousRelative. A mesh-owned AUTO_WRITE dictionary
uniform/hisaCFLStagnationState saves integral and previousRelative with the
flow checkpoint. History is optional: if absent at startup, integral starts
at zero and the first PTC call seeds previousRelative. Reading occurs at the
selected restart time before advancing Time, not lazily at the next step. NK does not update it; the next PTC call includes any intervening
NK residual change. Exponential attenuation is limited to exp(-700).
No extra residual evaluation or MPI collective. Gain must be finite >=0;
decay finite >0. Set gain0 to retain the existing damped SER behavior.
A fresh case must not copy controller history unless deliberately restarting.
These controls do not guarantee stagnation removal; a too-small CFL may also
cause slow convergence. Decay/gain defaults are not universal optimal settings.

Optional hybrid PTC trend gain scheduling (default disabled):
```
pseudoCoGainScheduling true;
pseudoCoTrendWeight    0.05;
pseudoCoTrendEpsilon   0.001;
pseudoCoMinGain        0.3;
pseudoCoStagnationGain 0;
pseudoCoLogDamping     1;
```
Use signed slope s=log(previousRelative)-log(relative), EMA m=(1-beta)*m+beta*s,
and q=qmin+(1-qmin)*max(m,0)/(max(m,0)+epsilon). Multiply the original SER target
by q before the existing rate and absolute limits. Epsilon is a user-set
positive transition scale (log residual change per PTC call), not an absolute
residual tolerance. At m=epsilon the correction is (1+qmin)/2; at m<=0 it is
qmin. It is a bounded schedule, not an accumulating integral or proven
convergence guarantee. Smaller epsilon releases the correction more readily.
Beta and qmin must be finite in (0,1], epsilon finite >0.
Cannot combine with a positive stagnation-integral gain. For pure scheduling
set log damping1 as shown. Non-hybrid native SER and NK remain unchanged.
Optional uniform/hisaCFLTrendState stores previousRelative, signed trend,
trendInitialized and diagnostic gain. Read at startup at the selected time;
without history the first call seeds previousRelative and uses m0=0, the
first available slope initializes the EMA. No extra residual or MPI collective.

### NK cubic line search

`nkLineSearch relaxedCubic` (default) follows the local ADflow `LSCubic`
acceptance policy: unit step, quadratic then cubic backtracking, Armijo
coefficient 0.01, and successive interpolation steps restricted to 0.1–0.5
of the previous step. `nkMinLineSearchStep 0.01` is the backtracking stopping
threshold, not a lower clamp. A finite admissible trial at or below that
threshold can be accepted without sufficient decrease. `strictCubic` retains
the former mandatory-decrease policy and 1e-8 backtracking cutoff.

A full-step SA residual norm increase above twice the base norm triggers
prelimiting (initial factor max(base/trial, 0.1), then halving). Invalid
trials also trigger preliminary halving, with at most ten preliminary
trials. The preliminary Armijo test follows ADflow's unscaled-slope RHS.
The SA norm is formed from the scaled SA residual entries; its fixed block
scale cancels in the growth ratio. This adds one scalar global reduction
per evaluated merit for the turbulence norm.

Intentional differences from ADflow: invalid/nonphysical trials are always
rejected; exhausting preliminary trials rejects and restores the base state;
we do not reset to a uniform field or change the user's NK switch threshold
on NaN. Existing linear-defect rejection and PTC fallback remain in force.
The `alpha` monitor includes `ls minimumStep` or `ls sufficientDecrease`,
and `prelimited` identifies preliminary backtracking. This is a primal
iteration policy; the residual and its adjoint definition are unchanged.

Standalone policy verification:
```
g++ -std=c++17 -I src/petscNK tests/petscNK/cubicSearchTest.C -o /tmp/cubicSearchTest
/tmp/cubicSearchTest
```

Only `relaxedCubic` and `strictCubic` are accepted line-search names.

### Stepwise residual feedback for PTC

Set `pseudoCoResidualFeedback` in `pseudoTime` to enable per-step feedback in
hybrid SER mode. The defaults are shown below. Disable the windowed plateau,
random, trend and stagnation controllers when using this option.

```foam
pseudoCoResidualFeedback
{
    enabled         true;
    targetDecay     0.002302585092994046;
    upThreshold     0.15;
    downThreshold   0.2302585092994046;
    increase        1.02;
    decrease        0.8;
    minFactor       1e-4;
    maxFactor       100;
    referenceWeight 0.1;
    trialTolerance  0.15;
    trialMaxSteps   40;
}
```

Each PTC update contributes its log residual decay minus `targetDecay` to
the growth integral, and the negated difference to the decrease integral.
Both integrals have a floor of zero. Crossing the growth threshold increases
the multiplier. Crossing the decrease threshold starts a trial if the new
multiplier would lower the actual CFL. If omitted, `downThreshold` defaults
to `100*targetDecay`.

The trial freezes the multiplier and compares the accumulated log residual
drop with the extrapolated pre-trial rate. `referenceWeight` is the weight of
the newest sample in that rate's exponential average. Improvement of at
least `trialTolerance` retains the multiplier. Deterioration by the same
tolerance or reaching `trialMaxSteps` without confirmation restores the
prior multiplier and waits for a change in the residual trend before another
decrease trial.

The multiplier applies to the uncapped SER target, followed by the existing
CFL rate and absolute bounds. With `pseudoCoLogDamping 1`, the rate limits
apply directly to the target-to-current CFL ratio. The multiplier, integrals
and pending trial are checkpointed in `uniform/hisaCFLFeedbackState`.
After an NK gap, the controller discards an unfinished trial and resets its
residual trend.
