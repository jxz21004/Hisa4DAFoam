#include "psiThermo.H"
#include "makeThermo.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "sensibleInternalEnergy.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "hePsiThermo.H"
#include "pureMixture.H"
#include "sutherlandConstPrTransport.H"
namespace Foam {
makeThermos(psiThermo, hePsiThermo, pureMixture, sutherlandConstPrTransport,
    sensibleInternalEnergy, hConstThermo, perfectGas, specie);
makeThermos(psiThermo, hePsiThermo, pureMixture, sutherlandConstPrTransport,
    sensibleEnthalpy, hConstThermo, perfectGas, specie);
}
