/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sutherlandConstPrTransport

Group
    grpSpecieTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    sutherlandConstPrTransportI.H
    sutherlandConstPrTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandConstPrTransport_H
#define sutherlandConstPrTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class sutherlandConstPrTransport;

template<class Thermo>
inline sutherlandConstPrTransport<Thermo> operator+
(
    const sutherlandConstPrTransport<Thermo>&,
    const sutherlandConstPrTransport<Thermo>&
);

template<class Thermo>
inline sutherlandConstPrTransport<Thermo> operator*
(
    const scalar,
    const sutherlandConstPrTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sutherlandConstPrTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class sutherlandConstPrTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sutherlandConstPrTransport
:
    public Thermo
{
    // Private Data

        // Sutherland's coefficients
        scalar As_, Ts_, Pr_ = 0.72;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from components
        inline sutherlandConstPrTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts,
            const scalar Pr = 0.72
        );

        //- Construct from two viscosities
        inline sutherlandConstPrTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline sutherlandConstPrTransport(const word&, const sutherlandConstPrTransport&);

        //- Construct from dictionary
        explicit sutherlandConstPrTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        sutherlandConstPrTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sutherlandConstPrTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<sutherlandConstPrTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherlandConstPr<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const sutherlandConstPrTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend sutherlandConstPrTransport operator+ <Thermo>
        (
            const sutherlandConstPrTransport&,
            const sutherlandConstPrTransport&
        );

        friend sutherlandConstPrTransport operator* <Thermo>
        (
            const scalar,
            const sutherlandConstPrTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sutherlandConstPrTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandConstPrTransportI.H"

#ifdef NoRepository
    #include "sutherlandConstPrTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
