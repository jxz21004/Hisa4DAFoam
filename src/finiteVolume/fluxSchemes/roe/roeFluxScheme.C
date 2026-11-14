// Former option name: roeCfl3dEntropyFix. Use roeAllWaveEntropyFix.
// GPL-3.0-or-later
#include "roeFluxScheme.H"
#include "roeKernel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "emptyFvPatch.H"
#include "zeroGradientFvPatchFields.H"
namespace Foam {
defineTypeNameAndDebug(roeFluxScheme,0);
addToRunTimeSelectionTable(fluxScheme,roeFluxScheme,dictionary);
roeFluxScheme::roeFluxScheme(const dictionary& dict,const fluidThermo& t,const volScalarField& rho,const volVectorField& U,const volVectorField&,const volScalarField& E)
:fluxScheme(typeName,dict,U.mesh()),mesh_(U.mesh()),thermo_(t),rho_(rho),U_(U),rhoE_(E),
 kappa_(dict.lookupOrDefault<scalar>("roeKappa",1.0/3.0)),entropy_(dict.lookupOrDefault<scalar>("roeEntropyFix",0.1)),limit_(dict.lookupOrDefault<Switch>("roeLimit",true)),allWaveEntropy_(dict.lookupOrDefault<Switch>("roeAllWaveEntropyFix",false)),
 limiter_(dict.lookupOrDefault<word>("roeLimiter",limit_ ? "hard" : "none")),
 smoothCoeff_(dict.lookupOrDefault<scalar>((limiter_=="vanAlbadaType" || limiter_=="vanAlbada") ? "roeVanAlbadaCoeff" : "roeVenkatCoeff",1)),
 smoothLength_(dict.lookupOrDefault<scalar>((limiter_=="vanAlbadaType" || limiter_=="vanAlbada") ? "roeVanAlbadaReferenceLength" : "roeVenkatReferenceLength",1)),
 smoothRho_(dict.lookupOrDefault<scalar>((limiter_=="vanAlbadaType" || limiter_=="vanAlbada") ? "roeVanAlbadaReferenceRho" : "roeVenkatReferenceRho",1)),
 smoothU_(dict.lookupOrDefault<scalar>((limiter_=="vanAlbadaType" || limiter_=="vanAlbada") ? "roeVanAlbadaReferenceU" : "roeVenkatReferenceU",1)),
 smoothP_(dict.lookupOrDefault<scalar>((limiter_=="vanAlbadaType" || limiter_=="vanAlbada") ? "roeVanAlbadaReferenceP" : "roeVenkatReferenceP",1)) {
    if(limiter_!="none" && limiter_!="hard" && limiter_!="venkatakrishnan" && limiter_!="vanAlbadaType" && limiter_!="vanAlbada")
        FatalErrorInFunction<<"roeLimiter must be none, hard, venkatakrishnan, vanAlbada or vanAlbadaType"<<exit(FatalError);
    if(dict.found("roeLimiter") && dict.found("roeLimit") && limit_!=(limiter_!="none"))
        FatalErrorInFunction<<"Conflicting roeLimit and roeLimiter settings"<<exit(FatalError);
    if((limiter_=="venkatakrishnan" || (limiter_=="vanAlbadaType" || limiter_=="vanAlbada")) &&
        (smoothCoeff_<=0 || smoothLength_<=0 || smoothRho_<=0 || smoothU_<=0 || smoothP_<=0))
        FatalErrorInFunction<<"Smooth limiter coefficient and reference scales must be positive"<<exit(FatalError);
    if(kappa_<0 || kappa_>1 || entropy_<0 || entropy_>1)
        FatalErrorInFunction<<"Roe requires 0<=roeKappa<=1 and 0<=roeEntropyFix<=1"<<exit(FatalError);
    if(mesh_.moving())FatalErrorInFunction<<"Roe supports stationary meshes only"<<exit(FatalError);
    forAll(mesh_.boundary(),p)if(mesh_.boundary()[p].coupled() && mesh_.boundary()[p].type() != "processor")FatalErrorInFunction<<"Roe supports processor interfaces, but not cyclic/AMI interfaces"<<exit(FatalError);
    Info<<"Roe UMUSCL: kappa="<<kappa_<<", entropy="<<entropy_<<", all-wave entropy="<<allWaveEntropy_<<", limiter="<<limiter_<<endl;
}
static scalar albadaRoe(const word& limiter,const scalar& delta,const scalar& a,
    const scalar& b,const scalar& eps2) {
    return limiter=="vanAlbada"
        ? hisaRoe::vanAlbadaIncrement(delta,a,b,eps2)
        : hisaRoe::vanAlbadaTypeIncrement(delta,a,b,eps2);
}
// Use the same facewise limiter for internal and processor faces.
static scalar limitRoe(const scalar& face,const scalar& cell,const scalar& low,
    const scalar& high,const word& limiter,const scalar& eps2) {
    if(limiter=="hard")return min(high,max(low,face));
    if(limiter=="none")return face;
    const scalar delta=face-cell;
    const scalar available=(delta>=0 ? high : low)-cell;
    return cell+hisaRoe::venkatIncrement(delta,available,eps2);
}
// Bound against the full adjacent-cell stencil, identically for both kappa values.
static void reconstructRoe(const volScalarField& q,const scalar k,const word& limiter,const volScalarField* smoothing,const scalar& reference,scalarField& l,scalarField& r,List<scalarField>& bl,List<scalarField>& br) {
    const fvMesh& m=q.mesh();auto grad=fvc::grad(q);
    grad.ref().correctBoundaryConditions();
    volScalarField low(IOobject("roeLowerBound",m.time().timeName(),m,IOobject::NO_READ,IOobject::NO_WRITE,false),q);
    volScalarField high(IOobject("roeUpperBound",m.time().timeName(),m,IOobject::NO_READ,IOobject::NO_WRITE,false),q);
    scalarField& lo=low.primitiveFieldRef();scalarField& hi=high.primitiveFieldRef();
    const auto& own=m.owner();const auto& nei=m.neighbour();
    const bool useBounds=limiter=="hard" || limiter=="venkatakrishnan";
    if(useBounds){
    forAll(own,f){label i=own[f],j=nei[f];lo[i]=min(lo[i],q[j]);hi[i]=max(hi[i],q[j]);lo[j]=min(lo[j],q[i]);hi[j]=max(hi[j],q[i]);}
    forAll(m.boundary(),p){const auto& bc=q.boundaryField()[p];const auto& cells=m.boundary()[p].faceCells();
        if(isA<emptyFvPatch>(m.boundary()[p]))continue;
        const scalarField neighbour(m.boundary()[p].coupled() ? bc.patchNeighbourField() : tmp<scalarField>(new scalarField(bc)));
        forAll(cells,f){label i=cells[f];lo[i]=min(lo[i],neighbour[f]);hi[i]=max(hi[i],neighbour[f]);}}
    // Exchange full-stencil bounds, not the patch interpolated primitive value.
    low.correctBoundaryConditions();high.correctBoundaryConditions();
    }
    l.setSize(own.size());r.setSize(own.size());
    forAll(own,f){label i=own[f],j=nei[f];vector d=m.C()[j]-m.C()[i];
        l[f]=hisaRoe::reconstruct(q[i],q[j],scalar(grad()[i]&(m.Cf()[f]-m.C()[i])),scalar(grad()[i]&d),k);
        r[f]=hisaRoe::reconstruct(q[j],q[i],scalar(grad()[j]&(m.Cf()[f]-m.C()[j])),scalar(-(grad()[j]&d)),k);
        if(limiter=="vanAlbadaType" || limiter=="vanAlbada"){
            const scalar jump=q[j]-q[i];
            const scalar ai=2*(grad()[i]&d)-jump;
            const scalar aj=-2*(grad()[j]&d)+jump;
            l[f]=q[i]+albadaRoe(limiter,scalar(l[f]-q[i]),ai,jump,scalar((*smoothing)[i]*sqr(reference)));
            r[f]=q[j]+albadaRoe(limiter,scalar(r[f]-q[j]),aj,scalar(-jump),scalar((*smoothing)[j]*sqr(reference)));
        } else if(limiter!="none"){
            const scalar ei=smoothing ? (*smoothing)[i]*sqr(reference) : scalar(0);
            const scalar ej=smoothing ? (*smoothing)[j]*sqr(reference) : scalar(0);
            l[f]=limitRoe(l[f],q[i],lo[i],hi[i],limiter,ei);
            r[f]=limitRoe(r[f],q[j],lo[j],hi[j],limiter,ej);}}
    bl.setSize(m.boundary().size());br.setSize(m.boundary().size());
    forAll(m.boundary(),p){const auto& patch=m.boundary()[p];if(!patch.coupled())continue;
        const auto& cells=patch.faceCells();
        const scalarField qn(q.boundaryField()[p].patchNeighbourField());
        const vectorField gn(grad().boundaryField()[p].patchNeighbourField());
        scalarField lon,hin;
        if(useBounds){lon=low.boundaryField()[p].patchNeighbourField();hin=high.boundaryField()[p].patchNeighbourField();}
        const vectorField delta(patch.delta());
        scalarField en;
        if(smoothing)en=smoothing->boundaryField()[p].patchNeighbourField();
        bl[p].setSize(patch.size());br[p].setSize(patch.size());
        forAll(cells,f){label i=cells[f];vector df=patch.Cf()[f]-m.C()[i];
            bl[p][f]=hisaRoe::reconstruct(q[i],qn[f],scalar(grad()[i]&df),scalar(grad()[i]&delta[f]),k);
            br[p][f]=hisaRoe::reconstruct(qn[f],q[i],scalar(gn[f]&(df-delta[f])),scalar(-(gn[f]&delta[f])),k);
            if(limiter=="vanAlbadaType" || limiter=="vanAlbada"){
                const scalar jump=qn[f]-q[i];
                const scalar ai=2*(grad()[i]&delta[f])-jump;
                const scalar aj=-2*(gn[f]&delta[f])+jump;
                bl[p][f]=q[i]+albadaRoe(limiter,scalar(bl[p][f]-q[i]),ai,jump,scalar((*smoothing)[i]*sqr(reference)));
                br[p][f]=qn[f]+albadaRoe(limiter,scalar(br[p][f]-qn[f]),aj,scalar(-jump),scalar(en[f]*sqr(reference)));
            } else if(limiter!="none"){
                const scalar ei=smoothing ? (*smoothing)[i]*sqr(reference) : scalar(0);
                const scalar ej=smoothing ? en[f]*sqr(reference) : scalar(0);
                bl[p][f]=limitRoe(bl[p][f],q[i],lo[i],hi[i],limiter,ei);
                br[p][f]=limitRoe(br[p][f],qn[f],lon[f],hin[f],limiter,ej);}
        }
    }
}
void roeFluxScheme::calcFlux(surfaceScalarField& phi,surfaceVectorField& phiUp,surfaceScalarField& phiEp,surfaceVectorField& Up) {
    const auto tg=thermo_.gamma();const auto& gamma=tg();const auto& p=thermo_.p();
    const scalar g=gMin(gamma.primitiveField());
    if(g<=1 || gMax(mag(gamma.primitiveField()-g))>1e-10)
        FatalErrorInFunction<<"Roe requires spatially constant gamma>1"<<exit(FatalError);
    // Reject unsupported caloric offsets rather than silently using a different energy equation.
    scalar error=0;
    forAll(rho_,i){scalar expected=p[i]/(g-1)+0.5*rho_[i]*magSqr(U_[i]);error=max(error,mag(rhoE_[i]-expected)/max(mag(expected),scalar(1)));}
    reduce(error,maxOp<scalar>());
    if(error>1e-7)FatalErrorInFunction<<"Roe requires calorically perfect gas with zero energy offset; mismatch="<<error<<exit(FatalError);
    // eps_q^2 = q_ref^2 (K h/L_ref)^3, h=cbrt(V). Keep geometry
    // active for shape AD; no global solution extrema or extra reductions.
    autoPtr<volScalarField> smoothing;
    if(limiter_=="venkatakrishnan" || (limiter_=="vanAlbadaType" || limiter_=="vanAlbada")){
        smoothing.reset(new volScalarField(
            IOobject("roeLimiterSmoothing",mesh_.time().timeName(),mesh_,IOobject::NO_READ,IOobject::NO_WRITE,false),
            mesh_,dimensionedScalar("zero",dimless,0),zeroGradientFvPatchScalarField::typeName));
        smoothing->primitiveFieldRef()=mesh_.V()*pow3(smoothCoeff_/smoothLength_);
        smoothing->correctBoundaryConditions();
    }
    const volScalarField* smoothPtr=smoothing.valid() ? &smoothing() : nullptr;
    scalarField l[5],r[5];List<scalarField> bl[5],br[5];
    reconstructRoe(rho_,kappa_,limiter_,smoothPtr,smoothRho_,l[0],r[0],bl[0],br[0]);reconstructRoe(p,kappa_,limiter_,smoothPtr,smoothP_,l[4],r[4],bl[4],br[4]);
    for(direction c=0;c<3;++c)reconstructRoe(U_.component(c)(),kappa_,limiter_,smoothPtr,smoothU_,l[c+1],r[c+1],bl[c+1],br[c+1]);
    label fallback=0;
    forAll(mesh_.owner(),f){label i=mesh_.owner()[f],j=mesh_.neighbour()[f];
        scalar ul[3],ur[3],n[3];vector nv=mesh_.Sf()[f]/mesh_.magSf()[f];
        for(int c=0;c<3;++c){ul[c]=l[c+1][f];ur[c]=r[c+1][f];n[c]=nv[c];}
        scalar rl=l[0][f],rr=r[0][f],pl=l[4][f],pr=r[4][f];
        if(rl<=0 || rr<=0 || pl<=0 || pr<=0){++fallback;rl=rho_[i];rr=rho_[j];pl=p[i];pr=p[j];for(int c=0;c<3;++c){ul[c]=U_[i][c];ur[c]=U_[j][c];}}
        if(rl<=0 || rr<=0 || pl<=0 || pr<=0)FatalErrorInFunction<<"Nonpositive cell state in Roe face "<<f<<exit(FatalError);
        auto fl=hisaRoe::flux(hisaRoe::conservative(rl,ul,pl,g),hisaRoe::conservative(rr,ur,pr,g),n,g,entropy_,allWaveEntropy_);
        const scalar area=mesh_.magSf()[f];phi[f]=area*fl.q[0];phiEp[f]=area*fl.q[4];
        for(int c=0;c<3;++c)phiUp[f][c]=area*fl.q[c+1];
    }
    // Physical boundary states already incorporate the characteristic/no-slip conditions.
    forAll(mesh_.boundary(),patch){if(isA<emptyFvPatch>(mesh_.boundary()[patch]))continue;
        if(mesh_.boundary()[patch].coupled()){
            const auto& cells=mesh_.boundary()[patch].faceCells();
            const scalarField rn(rho_.boundaryField()[patch].patchNeighbourField()),pn(p.boundaryField()[patch].patchNeighbourField());
            const vectorField un(U_.boundaryField()[patch].patchNeighbourField());
            forAll(cells,f){scalar ul[3],ur[3],n[3];label i=cells[f];
                scalar rl=bl[0][patch][f],rr=br[0][patch][f],pl=bl[4][patch][f],pr=br[4][patch][f];
                const scalar area=mesh_.magSf().boundaryField()[patch][f];vector nv=mesh_.Sf().boundaryField()[patch][f]/area;
                for(int c=0;c<3;++c){ul[c]=bl[c+1][patch][f];ur[c]=br[c+1][patch][f];n[c]=nv[c];}
                if(rl<=0 || rr<=0 || pl<=0 || pr<=0){++fallback;rl=rho_[i];rr=rn[f];pl=p[i];pr=pn[f];for(int c=0;c<3;++c){ul[c]=U_[i][c];ur[c]=un[f][c];}}
                if(rl<=0 || rr<=0 || pl<=0 || pr<=0)FatalErrorInFunction<<"Nonpositive Roe processor state"<<exit(FatalError);
                auto fl=hisaRoe::flux(hisaRoe::conservative(rl,ul,pl,g),hisaRoe::conservative(rr,ur,pr,g),n,g,entropy_,allWaveEntropy_);
                phi.boundaryFieldRef()[patch][f]=area*fl.q[0];phiEp.boundaryFieldRef()[patch][f]=area*fl.q[4];
                for(int c=0;c<3;++c)phiUp.boundaryFieldRef()[patch][f][c]=area*fl.q[c+1];
            }
            continue;
        }
        const auto& rb=rho_.boundaryField()[patch];const auto& ub=U_.boundaryField()[patch];const auto& pb=p.boundaryField()[patch];
        forAll(rb,f){scalar u[3],n[3],area=mesh_.magSf().boundaryField()[patch][f];vector nv=mesh_.Sf().boundaryField()[patch][f]/area;
            if(rb[f]<=0 || pb[f]<=0)FatalErrorInFunction<<"Nonpositive Roe boundary state"<<exit(FatalError);
            for(int c=0;c<3;++c){u[c]=ub[f][c];n[c]=nv[c];}
            auto fl=hisaRoe::physical(hisaRoe::conservative(rb[f],u,pb[f],g),n,g);
            phi.boundaryFieldRef()[patch][f]=area*fl.q[0];phiEp.boundaryFieldRef()[patch][f]=area*fl.q[4];
            for(int c=0;c<3;++c)phiUp.boundaryFieldRef()[patch][f][c]=area*fl.q[c+1];}}
    // Continuous viscous-work velocity, independent of the upwind/contact switch.
    Up=fvc::interpolate(U_)*mesh_.magSf();
    reduce(fallback,sumOp<label>());
    if(fallback)Info<<"Roe reconstruction positivity fallback face incidences (processor faces counted twice): "<<fallback<<endl;
}
}
