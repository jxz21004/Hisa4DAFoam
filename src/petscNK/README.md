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

PTC uses native HISA flow iteration and DAFoam SA correction. Its CFL continues
to follow the existing PTC configuration. NK does not use that CFL as a step
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
