# Roe entropy kernel regression

From the repository root:

```sh
g++ -std=c++17 -O2 -I src/finiteVolume/fluxSchemes/roe tests/roeEntropy/test_kernel.cpp -o /tmp/test_roe_entropy
/tmp/test_roe_entropy
```

The optional `roeCfl3dEntropyFix true` in fvSchemes applies Harten/Gnoffo
regularization to all wave families with epsilon = roeEntropyFix*(a+|u|+|v|+|w|).
Default false preserves the existing acoustic-only correction. The same option
is used in the frozen approximate Jacobian, without extra LHS amplification.
The test covers stationary-contact dissipation, conservation under normal
reversal, equal-state consistency, zero-coefficient equivalence and the C1 join.
It does not validate coupled solver convergence or adjoint gradients.
