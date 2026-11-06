# highOrderNewton regression support

Reusable HiSA/DAFoam checks: audit_connectivity.py, check_parity.py and check_derivatives.py. geometric_jacobian.py supplies the connectivity labels helper. dafoam/tools/validate_hisa_local.py calls these checks with explicit external case/library paths. Newton correction experiments and saved reports were archived. These checks validate frozen-state connectivity/derivatives, not converged aerodynamic gradients.

Project CFD validation belongs in `../../../validation/hisa_adflow_comparison`; large runtime cases belong under workspace `reports`.
