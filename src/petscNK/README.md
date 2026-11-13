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
Other monitor names: `ksp`, `linear` (recomputed linear defect), `alpha`.

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
