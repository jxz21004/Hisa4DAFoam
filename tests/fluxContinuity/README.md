# fluxContinuity regression support

prepare_probe.py, probe_state.py and verify_fix.py probe the energy residual across a face-flow-direction switch. The viscous work correction must remain continuous as mass flux passes through zero. Supply external case/reference/output paths; the shared wingAero runner is retained. Old numerical reports were archived.

Project CFD validation belongs in `../../../validation/hisa_adflow_comparison`; large runtime cases belong under workspace `reports`.
