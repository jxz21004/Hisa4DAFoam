#include "HisaPetscNK.H"
#include <cassert>
#include <iostream>
int main(int argc, char** argv)
{
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    {
        hisaNK::Array z(1, 1.), observed=z;
        auto residual=[&](const hisaNK::Array& x) {
            observed=x;
            return hisaNK::Array{ x[0]*x[0]-2. };
        };
        auto matrix=[&]() {
            Mat p; MatCreateAIJ(PETSC_COMM_WORLD,1,1,PETSC_DECIDE,PETSC_DECIDE,1,nullptr,0,nullptr,&p);
            PetscInt lo,hi; MatGetOwnershipRange(p,&lo,&hi);
            MatSetValue(p,lo,lo,2*z[0],INSERT_VALUES);
            MatAssemblyBegin(p,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(p,MAT_FINAL_ASSEMBLY);
            return p;
        };
        hisaNK::Options options; options.lag=1; options.ew=false;
        hisaNK::Solver solver(options,{1.},residual,[](const hisaNK::Array&){return true;},matrix);
        for(int i=0;i<4;++i) assert(solver.step(z).accepted);
        assert(std::abs(z[0]*z[0]-2.)<1e-9);
        z[0]=1.;
        hisaNK::Solver rejected(options,{1.},residual,[](const hisaNK::Array&){return false;},matrix);
        assert(!rejected.step(z).accepted);
        assert(z[0]==1. && observed[0]==1.);
        int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        if(!rank) std::cout<<"NK core: convergence and rejected-state restoration passed\n";
    }
    PetscFinalize();
}
