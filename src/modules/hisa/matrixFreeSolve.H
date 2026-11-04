// Optional steady, frozen-turbulence matrix-free flow solve. Off by default.
const dictionary& mfDict = dict.subOrEmptyDict(solverType);
if (!mfDict.lookupOrDefault<bool>("matrixFree", false))
{
    sol->solve(scalarVars_, vectorVars_, scalarResiduals_, vectorResiduals_, res_);
}
else
{
    if (!steadyState_ || mesh.moving() || solverType != "GMRES")
    {
        FatalErrorInFunction << "matrixFree requires steady GMRES on a fixed mesh"
            << exit(FatalError);
    }
    const scalar epsilon = mfDict.lookupOrDefault<scalar>("matrixFreeEpsilon", 1e-6);
    if (epsilon <= 0) FatalErrorInFunction << "matrixFreeEpsilon must be positive" << exit(FatalError);
    const scalarField frozenDdt(ddtCoeff);
    const volScalarField baseRho("mfBaseRho", rho), baseE("mfBaseRhoE", rhoE);
    const volVectorField baseMom("mfBaseRhoU", rhoU), baseU("mfBaseU", U);
    const volScalarField baseP("mfBaseP", p), baseHe("mfBaseHe", e), baseT("mfBaseT", T);
    const volScalarField baseR("mfBaseR", rhoR), baseER("mfBaseER", rhoER);
    const volVectorField baseMR("mfBaseMR", rhoUR);
    const scalar rs = max(gAverage(mag(rho.primitiveField())), scalar(1e-8));
    const scalar es = max(gAverage(mag(rhoE.primitiveField())), scalar(1));
    const scalar ms = max(gAverage(mag(rhoU.primitiveField())), rs*100);
    label evaluations = 0;
    bool checkedProduct = false;
    auto evaluate = [&](const scalar alpha, const scalarField& dr,
                        const vectorField& dm, const scalarField& de)
    {
        rho = baseRho; rhoU = baseMom; rhoE = baseE;
        p = baseP; e = baseHe; T = baseT; U = baseU;
        rho.primitiveFieldRef() += alpha*dr;
        rhoU.primitiveFieldRef() += alpha*dm;
        rhoE.primitiveFieldRef() += alpha*de;
        bounded_ = false;
        {
            #include "updateFields.H"
        }
        flux.calcFlux(phi, phiUp, phiEp, Up);
        {
            #include "residualsUpdate.H"
        }
        ddtCoeff = frozenDdt;
        ++evaluations;
    };
    auto merit = [&]() -> scalar
    {
        return sqrt(gSum(magSqr(rhoR.primitiveField()))/sqr(rs)
                  + gSum(magSqr(rhoUR.primitiveField()))/sqr(ms)
                  + gSum(magSqr(rhoER.primitiveField()))/sqr(es));
    };
    const scalar initialMerit = merit();
    const scalarField zeroS(mesh.nCells(), 0);
    const vectorField zeroV(mesh.nCells(), vector::zero);
    evaluate(0, zeroS, zeroV, zeroS);
    const scalar repeatError = sqrt
    (
        gSum(magSqr(rhoR.primitiveField()-baseR.primitiveField()))/sqr(rs)
      + gSum(magSqr(rhoUR.primitiveField()-baseMR.primitiveField()))/sqr(ms)
      + gSum(magSqr(rhoER.primitiveField()-baseER.primitiveField()))/sqr(es)
    )/max(initialMerit, scalar(VSMALL));
    if (repeatError > 1e-8)
    {
        FatalErrorInFunction << "matrixFree zero-step residual drift " << repeatError << exit(FatalError);
    }
    rhoR = baseR; rhoUR = baseMR; rhoER = baseER;
    sol->setMatrixProduct
    (
        [&](PtrList<volScalarField>& sv, PtrList<volVectorField>& vv,
            PtrList<volScalarField>& sr, PtrList<volVectorField>& vr)
        {
            const scalarField dr(sv[0].primitiveField()), de(sv[1].primitiveField());
            const vectorField dm(vv[0].primitiveField());
            const scalar vn = max(max(gMax(mag(dr))/rs, gMax(mag(de))/es), gMax(mag(dm))/ms);
            if (vn < VSMALL)
            {
                sr[0].primitiveFieldRef() = scalar(0); sr[1].primitiveFieldRef() = scalar(0);
                vr[0].primitiveFieldRef() = vector::zero;
                return;
            }
            const scalar h = epsilon/vn;
            evaluate(h, dr, dm, de);
            const scalarField plusR(rhoR.primitiveField()), plusE(rhoER.primitiveField());
            const vectorField plusM(rhoUR.primitiveField());
            evaluate(-h, dr, dm, de);
            // Native residual is -div(flux): A = pseudo-time diagonal - dR/dW.
            // matrixMul and the RHS are both in strong form (no cell volume factor).
            sr[0].primitiveFieldRef() = frozenDdt*dr - (plusR-rhoR.primitiveField())/(2*h);
            sr[1].primitiveFieldRef() = frozenDdt*de - (plusE-rhoER.primitiveField())/(2*h);
            vr[0].primitiveFieldRef() = frozenDdt*dm - (plusM-rhoUR.primitiveField())/(2*h);
            if (!checkedProduct && mfDict.lookupOrDefault<bool>("matrixFreeCheck", false))
            {
                const scalarField fullR((plusR-rhoR.primitiveField())/(2*h));
                const scalarField fullE((plusE-rhoER.primitiveField())/(2*h));
                const vectorField fullM((plusM-rhoUR.primitiveField())/(2*h));
                evaluate(0.5*h, dr, dm, de);
                const scalarField halfPlusR(rhoR.primitiveField()), halfPlusE(rhoER.primitiveField());
                const vectorField halfPlusM(rhoUR.primitiveField());
                evaluate(-0.5*h, dr, dm, de);
                const scalarField halfR((halfPlusR-rhoR.primitiveField())/h);
                const scalarField halfE((halfPlusE-rhoER.primitiveField())/h);
                const vectorField halfM((halfPlusM-rhoUR.primitiveField())/h);
                scalar err = sqrt(gSum(magSqr(fullR-halfR))/sqr(rs)
                    +gSum(magSqr(fullM-halfM))/sqr(ms)+gSum(magSqr(fullE-halfE))/sqr(es));
                scalar ref = sqrt(gSum(magSqr(halfR))/sqr(rs)
                    +gSum(magSqr(halfM))/sqr(ms)+gSum(magSqr(halfE))/sqr(es));
                Info << "MatrixFreeCheck: epsilon " << epsilon << " relativeDifference "
                     << err/max(ref,scalar(VSMALL)) << endl;
                checkedProduct = true;
            }
            evaluate(0, zeroS, zeroV, zeroS);
            rhoR = baseR; rhoUR = baseMR; rhoER = baseER;
        }
    );
    sol->solve(scalarVars_, vectorVars_, scalarResiduals_, vectorResiduals_, res_);
    const scalarField dr(rho.primitiveField()-baseRho.primitiveField());
    const scalarField de(rhoE.primitiveField()-baseE.primitiveField());
    const vectorField dm(rhoU.primitiveField()-baseMom.primitiveField());
    scalar alpha = 1;
    bool accepted = false;
    scalar finalMerit = initialMerit;
    const label maxTrials = mfDict.lookupOrDefault<label>("matrixFreeMaxTrials", 12);
    for (label trial=0; trial<maxTrials; ++trial)
    {
        scalarField trialRho(baseRho.primitiveField()+alpha*dr);
        vectorField trialMom(baseMom.primitiveField()+alpha*dm);
        scalarField trialE(baseE.primitiveField()+alpha*de);
        if (gMin(trialRho)>1e-3 && gMin(trialE - 0.5*magSqr(trialMom)/max(trialRho,scalar(1e-30)))>0)
        {
            evaluate(alpha, dr, dm, de);
            finalMerit = merit();
            if (!bounded_ && finalMerit <= (1-1e-4*alpha)*initialMerit)
            {
                accepted = true;
                break;
            }
        }
        alpha *= 0.5;
    }
    if (!accepted)
    {
        evaluate(0, zeroS, zeroV, zeroS);
        finalMerit = initialMerit;
        alpha = 0;
    }
    Info << "MatrixFreeFlow: alpha " << alpha << " meritRatio "
        << finalMerit/max(initialMerit,scalar(VSMALL)) << " residualEvaluations "
        << evaluations << " repeatError " << repeatError << endl;
}
