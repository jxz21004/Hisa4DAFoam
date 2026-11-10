# Fixed-reference flow residual convergence

Opt in through `pseudoTime` in fvSolution:

```
fixedInitialResidual true;
fixedResidualRelTol 1e-8;
fixedResidualAbsTol (0 0 0 0 0);
```

This replaces the outer GMRES-normalized convergence test for steady flow only.
Inner GMRES tolerances, CFL control and the underlying discretization are unchanged.
Global cell-count RMS is computed from the current strong-form physical residual,
in order rho, rhoUx, rhoUy, rhoUz, rhoE, before the next flow update.
Each component must satisfy either its absolute tolerance (in that component's
residual units) or the fixed initial relative tolerance. Set physically justified
absolute tolerances for near-zero components; an exactly-zero baseline accepts
only zero residual unless an absolute tolerance is supplied.

The log prints `FixedFlowResidualRMS`, `FixedFlowResidualRelative`, and
`FixedFlowResidualInitial` every iteration. These are pre-update residuals, unlike
function-object samples stamped after an update. No extra flux evaluation is added.

The initial RMS vector is stored as `fixedFlowInitialRMS` in the checkpoint's
`uniform/pseudotimeState` (in each processor directory for MPI). A restart must
retain this file. Restarting an older case without the baseline is an error unless
`allowNewResidualBaseline true` explicitly accepts a NEW reference at that restart.
This flag does not reset an existing reference. Never compare a new reference's
relative values as if they were relative to the original uniform flow.

The standalone executable must be rebuilt as well as libHisaSolver: when flow
convergence is detected it skips both the flow update and the subsequent turbulence
correction, and ends at the next control-loop check.

Scope: FIVE FLOW EQUATIONS ONLY. This does not calculate the unnormalized SA
steady-equation residual. The message explicitly says SA is not included; it is
not certification of full RANS/adjoint convergence. AD builds are not validated.
The legacy mode remains default; its latched convergence flag now rechecks the
latest solver residual instead of accepting a rebound unconditionally.
