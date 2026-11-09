// GPL-3.0-or-later
#ifndef roeFluxScheme_H
#define roeFluxScheme_H
#include "fluxScheme.H"
#include "fluidThermo.H"
namespace Foam {
class roeFluxScheme : public fluxScheme {
    const fvMesh& mesh_;
    const fluidThermo& thermo_;
    const volScalarField& rho_;
    const volVectorField& U_;
    const volScalarField& rhoE_;
    scalar kappa_, entropy_;
    bool limit_, cfl3dEntropy_;
public:
    TypeName("Roe");
    roeFluxScheme(const dictionary&,const fluidThermo&,const volScalarField&,const volVectorField&,const volVectorField&,const volScalarField&);
    virtual void calcFlux(surfaceScalarField&,surfaceVectorField&,surfaceScalarField&,surfaceVectorField&);
};
}
#endif
