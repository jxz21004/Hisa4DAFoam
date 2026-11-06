# Repository scope

Keep this repository focused on solver source, build tooling, upstream examples and reusable code regression tests.

- Project CFD comparisons, parameter sweeps, convergence studies and aerodynamic/adjoint experiments belong in `/mdouser/test_aerostructure/validation/hisa_adflow_comparison`.
- Runtime cases and large outputs belong under `/mdouser/test_aerostructure/reports`, never in this repository's `tests` or source directories.
- Retain reusable regression helpers needed by build/validation tools; use explicit external case/output paths.
- Preserve solver changes and installed build dependencies during cleanup. Do not revert scientific changes to make Git clean.
- Do not commit caches, generated logs, flow fields or archived experiment reports here.
