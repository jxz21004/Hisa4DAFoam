// GPL-3.0-or-later
#ifndef roeJacobian_H
#define roeJacobian_H
#include "laxFriedrichs.H"
namespace Foam {
class roeJacobian : public laxFriedrichs {
public:
    TypeName("Roe");
    roeJacobian(const dictionary&,const fvMesh&,const volScalarField&,const volVectorField&,const volScalarField&,const fluidThermo&);
    virtual void addInviscidJacobian(compressibleJacobianMatrix&);
};
}
#endif
