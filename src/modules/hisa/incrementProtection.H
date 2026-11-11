// Protect one conservative update, without changing the residual or its Jacobian.
// Pressure ratio below uses perfect-gas internal-energy density; reject other EOS.
const bool incrementProtection = solnControl_->dict().lookupOrDefault<bool>("incrementProtection", false);
if (incrementProtection)
{
    const dictionary flow=mesh.solutionDict().subDict("flowSolver").subOrEmptyDict("GMRES");
    if (flow.lookupOrDefault<bool>("matrixFree",false) || flow.lookupOrDefault<bool>("limitedUpdate",false))
        FatalErrorInFunction << "incrementProtection cannot be combined with matrixFree/limitedUpdate" << exit(FatalError);
    const label interval=solnControl_->dict().lookupOrDefault<label>("incrementProtectionReportInterval",100);
    const bool report=interval>0 && runTime.timeIndex()%interval==0;
    // Roe already validates the calorically perfect, zero-offset energy model.
    if (flux.type() != "Roe")
        FatalErrorInFunction << "incrementProtection currently validated for Roe only" << exit(FatalError);
    const auto& r0=scalarVarsPrevIter_[0];
    const auto& E0=scalarVarsPrevIter_[1];
    const auto& m0=vectorVarsPrevIter_[0];
    scalar minimumAlpha=1;
    label limited=0,backtracks=0,firstCell=-1;
    forAll(rho,i)
    {
        const scalar internal0=E0[i]-0.5*magSqr(m0[i])/r0[i];
        if (!(r0[i]>0 && internal0>0))
            FatalErrorInFunction << "Invalid base state in incrementProtection cell " << i << exit(FatalError);
        const scalar dr=rho[i]-r0[i],dE=rhoE[i]-E0[i];
        const vector dm=rhoU[i]-m0[i];
        if (!std::isfinite(dr) || !std::isfinite(dE) ||
            !std::isfinite(dm.x()) || !std::isfinite(dm.y()) || !std::isfinite(dm.z()))
            FatalErrorInFunction << "Nonfinite update in incrementProtection cell " << i << exit(FatalError);
        // Linear pressure change; final nonlinear pressure is checked below.
        const vector u0=m0[i]/r0[i];
        const scalar dpRel=(dE-(u0&dm)+0.5*magSqr(u0)*dr)/internal0;
        const scalar ar=hisaIncrementProtection::coefficient<scalar>(dr/r0[i]);
        const scalar ap=hisaIncrementProtection::coefficient<scalar>(dpRel);
        scalar alpha=min(ar,ap);
        if (ar<1 || ap<1) ++limited;
        bool accepted=false;
        for (label trial=0;trial<80;++trial)
        {
            const scalar rt=r0[i]+alpha*dr;
            const vector mt=m0[i]+alpha*dm;
            const scalar Et=E0[i]+alpha*dE;
            // Also cap actual nonlinear density/pressure loss to 50%.
            if (std::isfinite(rt) && std::isfinite(Et) && rt>=0.5*r0[i] &&
                Et-0.5*magSqr(mt)/rt>=0.5*internal0)
            {
                rho[i]=rt;rhoU[i]=mt;rhoE[i]=Et;accepted=true;break;
            }
            alpha*=0.5;++backtracks;
        }
        if (!accepted || !(alpha>0))
            FatalErrorInFunction << "incrementProtection cannot find admissible update in cell " << i << exit(FatalError);
        if (alpha<minimumAlpha){minimumAlpha=alpha;firstCell=i;}
    }
    if (report) {
        if (Pstream::parRun())
            Pout << "IncrementProtection local: minAlpha=" << minimumAlpha << " cell=" << firstCell << endl;
        reduce(minimumAlpha,minOp<scalar>());reduce(limited,sumOp<label>());reduce(backtracks,sumOp<label>());
        Info << "IncrementProtection: minAlpha=" << minimumAlpha << " limitedCells=" << limited
             << " backtracks=" << backtracks << " localMinCell=" << firstCell << endl;
    }
}
