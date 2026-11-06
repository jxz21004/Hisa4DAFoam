// GPL-3.0-or-later
#include "roeJacobian.H"
#include "roeKernel.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam {
defineTypeNameAndDebug(roeJacobian,0);
addToRunTimeSelectionTable(jacobianInviscid,roeJacobian,dictionary);
roeJacobian::roeJacobian(const dictionary& d,const fvMesh& m,const volScalarField& r,const volVectorField& u,const volScalarField& e,const fluidThermo& t)
:laxFriedrichs(d,m,r,u,e,t){}
template<class T> static void initializeRoe(fvjMatrix<T>& a,const fvMesh& m) {
    a.diag()=pTraits<T>::zero;a.upper()=pTraits<T>::zero;a.lower()=pTraits<T>::zero;
    a.interfacesUpper().resize(m.boundary().size());a.interfacesLower().resize(m.boundary().size());
    forAll(m.boundary(),p){a.interfacesUpper().set(p,new Field<T>(m.boundary()[p].size(),pTraits<T>::zero));a.interfacesLower().set(p,new Field<T>(m.boundary()[p].size(),pTraits<T>::zero));}
}
template<class T> static void insertRoe(fvjMatrix<T>& a,label f,label i,label j,const T& L,const T& R) {
    a.diag()[i]+=L;a.upper()[f]+=R;a.lower()[f]-=L;a.diag()[j]-=R;
}
template<class T> static void insertRoeInterface(fvjMatrix<T>& a,label patch,label f,label i,const T& L,const T& R) {
    a.diag()[i]+=L;
    a.interfacesUpper()[patch][f]=R;
    a.interfacesLower()[patch][f]=-L;
}
void roeJacobian::addInviscidJacobian(compressibleJacobianMatrix& a) {
    if(mesh_.moving())FatalErrorInFunction<<"Roe Jacobian supports stationary meshes only"<<exit(FatalError);
    forAll(mesh_.boundary(),p)if(mesh_.boundary()[p].coupled() && mesh_.boundary()[p].type() != "processor")FatalErrorInFunction<<"Roe Jacobian supports processor interfaces, but not cyclic/AMI interfaces"<<exit(FatalError);
    gamma_=thermo_.gamma();meshPhi_.clear();
    const auto& schemes=mesh_.thisDb().lookupObject<IOdictionary>("fvSchemes");
    scalar entropy=schemes.lookupOrDefault<scalar>("roeEntropyFix",0.1);
    auto& cR=a.dSByS(0,0);auto& cM=a.dSByV(0,0);auto& cE=a.dSByS(0,1);
    auto& mR=a.dVByS(0,0);auto& mM=a.dVByV(0,0);auto& mE=a.dVByS(0,1);
    auto& eR=a.dSByS(1,0);auto& eM=a.dSByV(1,0);auto& eE=a.dSByS(1,1);
    initializeRoe(cR,mesh_);initializeRoe(cM,mesh_);initializeRoe(cE,mesh_);
    initializeRoe(mR,mesh_);initializeRoe(mM,mesh_);initializeRoe(mE,mesh_);
    initializeRoe(eR,mesh_);initializeRoe(eM,mesh_);initializeRoe(eE,mesh_);
    forAll(mesh_.owner(),f){label i=mesh_.owner()[f],j=mesh_.neighbour()[f];scalar area=mesh_.magSf()[f];vector nv=mesh_.Sf()[f]/area;scalar n[3];
        hisaRoe::State<scalar> l,r;l.q[0]=rho_[i];r.q[0]=rho_[j];l.q[4]=rhoE_[i];r.q[4]=rhoE_[j];
        for(int c=0;c<3;++c){l.q[c+1]=rhoU_[i][c];r.q[c+1]=rhoU_[j][c];n[c]=nv[c];}
        scalar g=gamma_()[i];
        if(l.q[0]<=0 || r.q[0]<=0 || hisaRoe::pressure(l,g)<=0 || hisaRoe::pressure(r,g)<=0)
            FatalErrorInFunction<<"Nonpositive state in Roe Jacobian"<<exit(FatalError);
        scalar L[5][5],R[5][5];hisaRoe::frozenBlocks(l,r,n,g,entropy,L,R);
        for(int x=0;x<5;++x)for(int y=0;y<5;++y){L[x][y]*=area;R[x][y]*=area;}
        insertRoe(cR,f,i,j,L[0][0],R[0][0]);insertRoe(cE,f,i,j,L[0][4],R[0][4]);
        insertRoe(eR,f,i,j,L[4][0],R[4][0]);insertRoe(eE,f,i,j,L[4][4],R[4][4]);
        insertRoe(cM,f,i,j,vector(L[0][1],L[0][2],L[0][3]),vector(R[0][1],R[0][2],R[0][3]));
        insertRoe(eM,f,i,j,vector(L[4][1],L[4][2],L[4][3]),vector(R[4][1],R[4][2],R[4][3]));
        insertRoe(mR,f,i,j,vector(L[1][0],L[2][0],L[3][0]),vector(R[1][0],R[2][0],R[3][0]));
        insertRoe(mE,f,i,j,vector(L[1][4],L[2][4],L[3][4]),vector(R[1][4],R[2][4],R[3][4]));
        insertRoe(mM,f,i,j,tensor(L[1][1],L[1][2],L[1][3],L[2][1],L[2][2],L[2][3],L[3][1],L[3][2],L[3][3]),tensor(R[1][1],R[1][2],R[1][3],R[2][1],R[2][2],R[2][3],R[3][1],R[3][2],R[3][3]));
    }
    // Each rank contributes its local row. The opposite rank supplies the other
    // row with reversed normal; interfaceUpper is +R in jacobianMatrix::matrixVectorMult.
    forAll(mesh_.boundary(),patch){const auto& bp=mesh_.boundary()[patch];if(!bp.coupled())continue;
        const scalarField rn(rho_.boundaryField()[patch].patchNeighbourField());
        const vectorField mn(rhoU_.boundaryField()[patch].patchNeighbourField());
        const scalarField en(rhoE_.boundaryField()[patch].patchNeighbourField());
        const auto& cells=bp.faceCells();
        forAll(cells,f){label i=cells[f];scalar area=mesh_.magSf().boundaryField()[patch][f];vector nv=mesh_.Sf().boundaryField()[patch][f]/area;scalar n[3];
            hisaRoe::State<scalar> l,r;l.q[0]=rho_[i];r.q[0]=rn[f];l.q[4]=rhoE_[i];r.q[4]=en[f];
            for(int c=0;c<3;++c){l.q[c+1]=rhoU_[i][c];r.q[c+1]=mn[f][c];n[c]=nv[c];}
            scalar g=gamma_()[i];
            if(l.q[0]<=0 || r.q[0]<=0 || hisaRoe::pressure(l,g)<=0 || hisaRoe::pressure(r,g)<=0)
                FatalErrorInFunction<<"Nonpositive processor state in Roe Jacobian"<<exit(FatalError);
            scalar L[5][5],R[5][5];hisaRoe::frozenBlocks(l,r,n,g,entropy,L,R);
            for(int x=0;x<5;++x)for(int y=0;y<5;++y){L[x][y]*=area;R[x][y]*=area;}
        insertRoeInterface(cR,patch,f,i,L[0][0],R[0][0]);insertRoeInterface(cE,patch,f,i,L[0][4],R[0][4]);
        insertRoeInterface(eR,patch,f,i,L[4][0],R[4][0]);insertRoeInterface(eE,patch,f,i,L[4][4],R[4][4]);
        insertRoeInterface(cM,patch,f,i,vector(L[0][1],L[0][2],L[0][3]),vector(R[0][1],R[0][2],R[0][3]));
        insertRoeInterface(eM,patch,f,i,vector(L[4][1],L[4][2],L[4][3]),vector(R[4][1],R[4][2],R[4][3]));
        insertRoeInterface(mR,patch,f,i,vector(L[1][0],L[2][0],L[3][0]),vector(R[1][0],R[2][0],R[3][0]));
        insertRoeInterface(mE,patch,f,i,vector(L[1][4],L[2][4],L[3][4]),vector(R[1][4],R[2][4],R[3][4]));
        insertRoeInterface(mM,patch,f,i,tensor(L[1][1],L[1][2],L[1][3],L[2][1],L[2][2],L[2][3],L[3][1],L[3][2],L[3][3]),tensor(R[1][1],R[1][2],R[1][3],R[2][1],R[2][2],R[2][3],R[3][1],R[3][2],R[3][3]));
        }
    }
    // Physical-boundary chain rules are inherited; all boundary fluxes use the BC state.
    // Reconstruction derivatives, Roe-average derivatives and SA coupling are intentionally frozen.
}
}
