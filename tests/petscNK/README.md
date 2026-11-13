# PETSc NK kernel regression

With MPI and PETSc environment loaded:

```sh
mpicxx -std=c++17 test_core.cpp -I../../src/petscNK \
  -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include \
  -L$PETSC_DIR/$PETSC_ARCH/lib -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib \
  -lpetsc -o /tmp/hisa_nk_core_test
mpirun -np 8 /tmp/hisa_nk_core_test
```

Checks a known nonlinear root and that an inadmissible update leaves both the
state vector and callback's physical state at the base point. CFD integration
experiments live outside this repository.

Also checks a constant nonzero residual (zero Jacobian): exactly two PC builds,
one retry, rejected status marked stagnant, and restoration of the base state.
A two-scale linear system also checks recovery: with a one-step Krylov cap,
the first identity PC stagnates, a refreshed exact diagonal PC succeeds, and
only the successful retry updates the state.
