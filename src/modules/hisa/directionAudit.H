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
const scalarField dr(rho.primitiveField()-bRho.primitiveField()), de(rhoE.primitiveField()-bE.primitiveField());
const vectorField dm(rhoU.primitiveField()-bMom.primitiveField());
PtrList<volScalarField> ds(2), ar(2);
PtrList<volVectorField> dv(1), am(1);
ds.set(0,new volScalarField("auditDr",bRho)); ds.set(1,new volScalarField("auditDe",bE));
dv.set(0,new volVectorField("auditDm",bMom));
ar.set(0,new volScalarField("auditAr",bR)); ar.set(1,new volScalarField("auditAe",bER));
am.set(0,new volVectorField("auditAm",bMR));
ds[0].primitiveFieldRef()=dr; ds[1].primitiveFieldRef()=de; dv[0].primitiveFieldRef()=dm;
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
         << " trueDerivativeNorm " << norm(tr,tm,te) << endl;
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
