# gmresEarlyExit regression support

run_regression.py checks early exit, restart-budget exhaustion and a single Krylov direction using the shared wingAero/run_case.py helper. Supply --template, --library-dir and --output outside this repository. This tests the linear solver stopping rule, not nonlinear CFD convergence.

Project CFD validation belongs in `../../../validation/hisa_adflow_comparison`; large runtime cases belong under workspace `reports`.
