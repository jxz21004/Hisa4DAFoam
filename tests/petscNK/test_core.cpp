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
        int builds=0;
        auto constantResidual=[&](const hisaNK::Array& x) {
            observed=x;
            return hisaNK::Array{1.};
        };
        auto countedMatrix=[&]() {++builds;return matrix();};
        hisaNK::Solver stagnant(options,{1.},constantResidual,
            [](const hisaNK::Array&){return true;},countedMatrix);
        auto stalled=stagnant.step(z);
        assert(!stalled.accepted && stalled.stagnated && stalled.retries==1);
        assert(builds==2 && z[0]==1. && observed[0]==1.);
        // A stale identity PC cannot solve this two-scale system in one step;
        // the rebuilt exact diagonal PC must recover without changing the base.
        hisaNK::Array linearZ{1., .01}; int linearBuilds=0;
        auto twoScale=[](const hisaNK::Array& x) {return hisaNK::Array{x[0],100*x[1]};};
        auto refreshingPC=[&]() {
            ++linearBuilds; Mat p;
            MatCreateAIJ(PETSC_COMM_WORLD,2,2,PETSC_DECIDE,PETSC_DECIDE,1,nullptr,0,nullptr,&p);
            PetscInt lo,hi;MatGetOwnershipRange(p,&lo,&hi);
            MatSetValue(p,lo,lo,1.,INSERT_VALUES);
            MatSetValue(p,lo+1,lo+1,linearBuilds==1?1.:100.,INSERT_VALUES);
            MatAssemblyBegin(p,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(p,MAT_FINAL_ASSEMBLY);
            return p;
        };
        hisaNK::Options limited=options;limited.maxLinear=1;limited.subspace=1;
        limited.stagnationTol=.5;
        hisaNK::Solver recover(limited,{1.,1.},twoScale,
            [](const hisaNK::Array&){return true;},refreshingPC);
        auto recovered=recover.step(linearZ);
        assert(recovered.accepted && recovered.retries==1 && !recovered.stagnated);
        assert(linearBuilds==2 && hisaNK::norm(twoScale(linearZ))<1e-6);
        int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        if(!rank) std::cout<<"NK core: convergence, restoration and bounded stagnation retry passed\n";
    }
    PetscFinalize();
}
