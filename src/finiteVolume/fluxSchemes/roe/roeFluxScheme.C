// GPL-3.0-or-later
#include "roeFluxScheme.H"
#include "roeKernel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "emptyFvPatch.H"
namespace Foam {
defineTypeNameAndDebug(roeFluxScheme,0);
addToRunTimeSelectionTable(fluxScheme,roeFluxScheme,dictionary);
roeFluxScheme::roeFluxScheme(const dictionary& dict,const fluidThermo& t,const volScalarField& rho,const volVectorField& U,const volVectorField&,const volScalarField& E)
:fluxScheme(typeName,dict,U.mesh()),mesh_(U.mesh()),thermo_(t),rho_(rho),U_(U),rhoE_(E),
 kappa_(dict.lookupOrDefault<scalar>("roeKappa",1.0/3.0)),entropy_(dict.lookupOrDefault<scalar>("roeEntropyFix",0.1)),limit_(dict.lookupOrDefault<Switch>("roeLimit",true)) {
    if(kappa_<0 || kappa_>1 || entropy_<0 || entropy_>1)
        FatalErrorInFunction<<"Roe requires 0<=roeKappa<=1 and 0<=roeEntropyFix<=1"<<exit(FatalError);
    if(Pstream::parRun() || mesh_.moving())FatalErrorInFunction<<"Experimental Roe supports serial stationary meshes only"<<exit(FatalError);
    forAll(mesh_.boundary(),p)if(mesh_.boundary()[p].coupled())FatalErrorInFunction<<"Experimental Roe does not yet support coupled patches"<<exit(FatalError);
    Info<<"Roe UMUSCL: kappa="<<kappa_<<", entropy="<<entropy_<<", bounded="<<limit_<<endl;
}
// Bound against the full adjacent-cell stencil, identically for both kappa values.
static void reconstructRoe(const volScalarField& q,const scalar k,const bool limit,scalarField& l,scalarField& r) {
    const fvMesh& m=q.mesh();const auto grad=fvc::grad(q);scalarField lo(q.primitiveField()),hi(lo);
    const auto& own=m.owner();const auto& nei=m.neighbour();
    forAll(own,f){label i=own[f],j=nei[f];lo[i]=min(lo[i],q[j]);hi[i]=max(hi[i],q[j]);lo[j]=min(lo[j],q[i]);hi[j]=max(hi[j],q[i]);}
    forAll(m.boundary(),p){const auto& bc=q.boundaryField()[p];const auto& cells=m.boundary()[p].faceCells();
        if(isA<emptyFvPatch>(m.boundary()[p]))continue;
        forAll(cells,f){label i=cells[f];lo[i]=min(lo[i],bc[f]);hi[i]=max(hi[i],bc[f]);}}
    l.setSize(own.size());r.setSize(own.size());
    forAll(own,f){label i=own[f],j=nei[f];vector d=m.C()[j]-m.C()[i];
        l[f]=hisaRoe::reconstruct(q[i],q[j],scalar(grad()[i]&(m.Cf()[f]-m.C()[i])),scalar(grad()[i]&d),k);
        r[f]=hisaRoe::reconstruct(q[j],q[i],scalar(grad()[j]&(m.Cf()[f]-m.C()[j])),scalar(-(grad()[j]&d)),k);
        if(limit){l[f]=min(hi[i],max(lo[i],l[f]));r[f]=min(hi[j],max(lo[j],r[f]));}}
}
void roeFluxScheme::calcFlux(surfaceScalarField& phi,surfaceVectorField& phiUp,surfaceScalarField& phiEp,surfaceVectorField& Up) {
    const auto tg=thermo_.gamma();const auto& gamma=tg();const auto& p=thermo_.p();
    const scalar g=gamma[0];
    if(g<=1 || gMax(mag(gamma.primitiveField()-g))>1e-10)
        FatalErrorInFunction<<"Roe requires spatially constant gamma>1"<<exit(FatalError);
    // Reject unsupported caloric offsets rather than silently using a different energy equation.
    scalar error=0;
    forAll(rho_,i){scalar expected=p[i]/(g-1)+0.5*rho_[i]*magSqr(U_[i]);error=max(error,mag(rhoE_[i]-expected)/max(mag(expected),scalar(1)));}
    if(error>1e-7)FatalErrorInFunction<<"Roe requires calorically perfect gas with zero energy offset; mismatch="<<error<<exit(FatalError);
    scalarField l[5],r[5];reconstructRoe(rho_,kappa_,limit_,l[0],r[0]);reconstructRoe(p,kappa_,limit_,l[4],r[4]);
    for(direction c=0;c<3;++c)reconstructRoe(U_.component(c)(),kappa_,limit_,l[c+1],r[c+1]);
    label fallback=0;
    forAll(mesh_.owner(),f){label i=mesh_.owner()[f],j=mesh_.neighbour()[f];
        scalar ul[3],ur[3],n[3];vector nv=mesh_.Sf()[f]/mesh_.magSf()[f];
        for(int c=0;c<3;++c){ul[c]=l[c+1][f];ur[c]=r[c+1][f];n[c]=nv[c];}
        scalar rl=l[0][f],rr=r[0][f],pl=l[4][f],pr=r[4][f];
        if(rl<=0 || rr<=0 || pl<=0 || pr<=0){++fallback;rl=rho_[i];rr=rho_[j];pl=p[i];pr=p[j];for(int c=0;c<3;++c){ul[c]=U_[i][c];ur[c]=U_[j][c];}}
        if(rl<=0 || rr<=0 || pl<=0 || pr<=0)FatalErrorInFunction<<"Nonpositive cell state in Roe face "<<f<<exit(FatalError);
        auto fl=hisaRoe::flux(hisaRoe::conservative(rl,ul,pl,g),hisaRoe::conservative(rr,ur,pr,g),n,g,entropy_);
        const scalar area=mesh_.magSf()[f];phi[f]=area*fl.q[0];phiEp[f]=area*fl.q[4];
        for(int c=0;c<3;++c)phiUp[f][c]=area*fl.q[c+1];
    }
    // Physical boundary states already incorporate the characteristic/no-slip conditions.
    forAll(mesh_.boundary(),patch){if(isA<emptyFvPatch>(mesh_.boundary()[patch]))continue;
        const auto& rb=rho_.boundaryField()[patch];const auto& ub=U_.boundaryField()[patch];const auto& pb=p.boundaryField()[patch];
        forAll(rb,f){scalar u[3],n[3],area=mesh_.magSf().boundaryField()[patch][f];vector nv=mesh_.Sf().boundaryField()[patch][f]/area;
            if(rb[f]<=0 || pb[f]<=0)FatalErrorInFunction<<"Nonpositive Roe boundary state"<<exit(FatalError);
            for(int c=0;c<3;++c){u[c]=ub[f][c];n[c]=nv[c];}
            auto fl=hisaRoe::physical(hisaRoe::conservative(rb[f],u,pb[f],g),n,g);
            phi.boundaryFieldRef()[patch][f]=area*fl.q[0];phiEp.boundaryFieldRef()[patch][f]=area*fl.q[4];
            for(int c=0;c<3;++c)phiUp.boundaryFieldRef()[patch][f][c]=area*fl.q[c+1];}}
    // Continuous viscous-work velocity, independent of the upwind/contact switch.
    Up=fvc::interpolate(U_)*mesh_.magSf();
    if(fallback)Info<<"Roe reconstruction positivity fallback faces: "<<fallback<<endl;
}
}
