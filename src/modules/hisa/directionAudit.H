// Diagnostic only: compare the assembled operator with true residual changes.
// Preserve the original flow update and stale flux so the iteration is unchanged.
if (!steadyState_ || mesh.moving() || solverType != "GMRES")
    FatalErrorInFunction << "directionAudit requires steady fixed-mesh GMRES" << exit(FatalError);
const volScalarField bRho("auditRho",rho), bE("auditE",rhoE), bP("auditP",p), bHe("auditHe",e), bT("auditT",T);
const volVectorField bMom("auditMom",rhoU), bU("auditU",U);
const volScalarField bR("auditR",rhoR), bER("auditER",rhoER);
const volVectorField bMR("auditMR",rhoUR);
const surfaceScalarField bPhi("auditPhi",phi), bPhiE("auditPhiE",phiEp);
const surfaceVectorField bPhiU("auditPhiU",phiUp), bUp("auditUp",Up);
const scalarField ddt(ddtCoeff);
const bool oldBounded = bounded_;
sol->solve(scalarVars_, vectorVars_, scalarResiduals_, vectorResiduals_, res_);
const volScalarField solvedRho("auditSolvedRho",rho), solvedE("auditSolvedE",rhoE);
const volVectorField solvedMom("auditSolvedMom",rhoU);
scalarField dr(rho.primitiveField()-bRho.primitiveField()), de(rhoE.primitiveField()-bE.primitiveField());
vectorField dm(rhoU.primitiveField()-bMom.primitiveField());
const bool storedDirection = mfDict.lookupOrDefault<bool>("directionAuditUseStored", false);
if (storedDirection)
{
    const volScalarField rd(IOobject("auditDeltaRho",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    const volScalarField ed(IOobject("auditDeltaE",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    const volVectorField md(IOobject("auditDeltaMom",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    dr=rd.primitiveField(); de=ed.primitiveField(); dm=md.primitiveField();
}
PtrList<volScalarField> ds(2), ar(2);
PtrList<volVectorField> dv(1), am(1);
ds.set(0,new volScalarField("auditDr",bRho)); ds.set(1,new volScalarField("auditDe",bE));
dv.set(0,new volVectorField("auditDm",bMom));
ar.set(0,new volScalarField("auditAr",bR)); ar.set(1,new volScalarField("auditAe",bER));
am.set(0,new volVectorField("auditAm",bMR));
ds[0].primitiveFieldRef()=dr; ds[1].primitiveFieldRef()=de; dv[0].primitiveFieldRef()=dm;
if (!storedDirection)
{
    ds[0].rename("auditDeltaRho"); ds[1].rename("auditDeltaE"); dv[0].rename("auditDeltaMom");
    ds[0].instance()=runTime.constant(); ds[1].instance()=runTime.constant(); dv[0].instance()=runTime.constant();
    ds[0].write(); ds[1].write(); dv[0].write();
}
sol->matrixProduct(ds,dv,ar,am);
const scalarField jr(ddt*dr-ar[0].primitiveField()), je(ddt*de-ar[1].primitiveField());
const vectorField jm(ddt*dm-am[0].primitiveField());
const scalar rs=max(gAverage(mag(bRho.primitiveField())),scalar(1e-8));
const scalar es=max(gAverage(mag(bE.primitiveField())),scalar(1));
const scalar ms=max(gAverage(mag(bMom.primitiveField())),scalar(rs*100));
auto norm = [&](const scalarField& r,const vectorField& m,const scalarField& en)->scalar
{ return sqrt(gSum(magSqr(r))/sqr(rs)+gSum(magSqr(m))/sqr(ms)+gSum(magSqr(en))/sqr(es)); };
const scalar initial=norm(bR.primitiveField(),bMR.primitiveField(),bER.primitiveField());
auto evaluate = [&](scalar alpha)
{
    rho=bRho; rhoE=bE; rhoU=bMom; p=bP; e=bHe; T=bT; U=bU;
    rho.primitiveFieldRef()+=alpha*dr; rhoE.primitiveFieldRef()+=alpha*de; rhoU.primitiveFieldRef()+=alpha*dm;
    bounded_=false;
    {
        #include "updateFields.H"
    }
    flux.calcFlux(phi,phiUp,phiEp,Up);
    {
        #include "residualsUpdate.H"
    }
    ddtCoeff=ddt;
};
evaluate(0);
Info << "DirectionAuditRepeat: " << norm(rhoR.primitiveField()-bR.primitiveField(),rhoUR.primitiveField()-bMR.primitiveField(),rhoER.primitiveField()-bER.primitiveField())/max(initial,scalar(VSMALL)) << endl;
for (scalar h : {scalar(1e-3),scalar(5e-4)})
{
    evaluate(h);
    const scalarField pr(rhoR.primitiveField()), pe(rhoER.primitiveField());
    const vectorField pm(rhoUR.primitiveField());
    evaluate(-h);
    const scalarField tr((pr-rhoR.primitiveField())/(2*h)), te((pe-rhoER.primitiveField())/(2*h));
    const vectorField tm((pm-rhoUR.primitiveField())/(2*h));
    Info << "DirectionAuditJ: h " << h << " relativeError " << norm(jr-tr,jm-tm,je-te)/max(norm(tr,tm,te),scalar(VSMALL))
         << " trueDerivativeNorm " << norm(tr,tm,te)
         << " densityError " << sqrt(gSum(magSqr(jr-tr)))/max(sqrt(gSum(magSqr(tr))),scalar(VSMALL))
         << " momentumError " << sqrt(gSum(magSqr(jm-tm)))/max(sqrt(gSum(magSqr(tm))),scalar(VSMALL))
         << " energyError " << sqrt(gSum(magSqr(je-te)))/max(sqrt(gSum(magSqr(te))),scalar(VSMALL)) << endl;
}
for (scalar alpha : {scalar(1),scalar(0.5),scalar(0.1),scalar(0.01)})
{
    evaluate(alpha);
    Info << "DirectionAuditStep: alpha " << alpha << " predictedRatio "
         << norm(bR.primitiveField()+alpha*jr,bMR.primitiveField()+alpha*jm,bER.primitiveField()+alpha*je)/max(initial,scalar(VSMALL))
         << " actualRatio " << norm(rhoR.primitiveField(),rhoUR.primitiveField(),rhoER.primitiveField())/max(initial,scalar(VSMALL)) << endl;
}
rho=solvedRho; rhoE=solvedE; rhoU=solvedMom; p=bP; e=bHe; T=bT; U=bU;
phi=bPhi; phiEp=bPhiE; phiUp=bPhiU; Up=bUp;
rhoR=bR; rhoER=bER; rhoUR=bMR; ddtCoeff=ddt; bounded_=oldBounded;
