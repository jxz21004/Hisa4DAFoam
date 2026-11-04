// Optional bounded or matrix-free steady flow solve. Both off by default.
const dictionary& mfDict = dict.subOrEmptyDict(solverType);
const bool matrixFree = mfDict.lookupOrDefault<bool>("matrixFree", false);
const bool limitedUpdate = mfDict.lookupOrDefault<bool>("limitedUpdate", false);
if (!matrixFree && !limitedUpdate)
{
    if (mfDict.lookupOrDefault<bool>("directionAudit", false))
    {
        #include "directionAudit.H"
    }
    else
    {
    sol->solve(scalarVars_, vectorVars_, scalarResiduals_, vectorResiduals_, res_);
    }
}
else
{
    if (!steadyState_ || mesh.moving() || solverType != "GMRES")
    {
        FatalErrorInFunction << "matrixFree/limitedUpdate requires steady GMRES on a fixed mesh"
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
    if (matrixFree)
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
    const scalar updateLimit = mfDict.lookupOrDefault<scalar>("updateLimit", 0.2);
    if (limitedUpdate)
    {
        if (updateLimit <= 0 || updateLimit >= 1)
            FatalErrorInFunction << "updateLimit must be between zero and one" << exit(FatalError);
        const scalar change = max
        (
            max(gMax(mag(dr)/max(mag(baseRho.primitiveField()),scalar(rs*1e-3))),
                gMax(mag(de)/max(mag(baseE.primitiveField()),scalar(es*1e-3)))),
            gMax(mag(dm)/max(mag(baseMom.primitiveField()),scalar(ms*0.1)))
        );
        alpha = min(scalar(1), updateLimit/max(change,scalar(VSMALL)));
    }
    bool accepted = false;
    scalar finalMerit = initialMerit;
    const label maxTrials = limitedUpdate
        ? mfDict.lookupOrDefault<label>("updateMaxTrials", 4)
        : mfDict.lookupOrDefault<label>("matrixFreeMaxTrials", 12);
    scalar pseudoMerit = initialMerit;
    const scalarField baseInternal(baseE.primitiveField()/baseRho.primitiveField()
        -0.5*magSqr(baseMom.primitiveField()/baseRho.primitiveField()));
    for (label trial=0; trial<maxTrials; ++trial)
    {
        scalarField trialRho(baseRho.primitiveField()+alpha*dr);
        vectorField trialMom(baseMom.primitiveField()+alpha*dm);
        scalarField trialE(baseE.primitiveField()+alpha*de);
        if (gMin(trialRho)>1e-3 && gMin(trialE - 0.5*magSqr(trialMom)/max(trialRho,scalar(1e-30)))>0)
        {
            const scalarField trialInternal(trialE/trialRho-0.5*magSqr(trialMom/trialRho));
            const bool physicalLimit = !limitedUpdate ||
                gMax(mag(trialInternal-baseInternal)/max(mag(baseInternal),scalar(1))) <= updateLimit;
            if (!physicalLimit) { alpha *= 0.5; continue; }
            evaluate(alpha, dr, dm, de);
            finalMerit = merit();
            pseudoMerit = sqrt
            (
                gSum(magSqr(alpha*frozenDdt*dr-rhoR.primitiveField()))/sqr(rs)
              + gSum(magSqr(alpha*frozenDdt*dm-rhoUR.primitiveField()))/sqr(ms)
              + gSum(magSqr(alpha*frozenDdt*de-rhoER.primitiveField()))/sqr(es)
            );
            const scalar acceptanceMerit = limitedUpdate ? pseudoMerit : finalMerit;
            if (!bounded_ && acceptanceMerit <= (1-1e-4*alpha)*initialMerit)
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
    if (limitedUpdate)
    {
        const scalar oldCfl = pseudoCoNum_().value();
        scalar nextCfl = oldCfl;
        if (!accepted || alpha < 0.25) nextCfl *= 0.5;
        else if (alpha >= 0.8 && finalMerit < 0.99*initialMerit) nextCfl *= 1.2;
        pseudoCoNum_().value() = max(pseudoCoNumMin_, min(pseudoCoNumMax_, nextCfl));
        Info << "LimitedUpdate: alpha " << alpha << " steadyRatio "
            << finalMerit/max(initialMerit,scalar(VSMALL)) << " pseudoRatio "
            << (accepted ? pseudoMerit/ max(initialMerit,scalar(VSMALL)) : scalar(1))
            << " residualEvaluations " << evaluations << " cfl " << oldCfl
            << " nextCfl " << pseudoCoNum_().value() << endl;
    }
    Info << (matrixFree ? "MatrixFreeFlow: alpha " : "LimitedUpdateCheck: alpha ") << alpha << " meritRatio "
        << finalMerit/max(initialMerit,scalar(VSMALL)) << " residualEvaluations "
        << evaluations << " repeatError " << repeatError << endl;
}
