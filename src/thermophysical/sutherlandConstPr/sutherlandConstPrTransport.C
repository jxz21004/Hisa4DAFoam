/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "sutherlandConstPrTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::sutherlandConstPrTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").get<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandConstPrTransport<Thermo>::sutherlandConstPrTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict)),
    Pr_(dict.subDict("transport").getOrDefault<scalar>("Pr", 0.72))
{
    if (!(Pr_ > 0)) FatalErrorInFunction << "Pr must be positive" << exit(FatalError);
}


template<class Thermo>
Foam::sutherlandConstPrTransport<Thermo>::sutherlandConstPrTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict)),
    Pr_(dict.subDict("transport").getOrDefault<scalar>("Pr", 0.72))
{
    if (!(Pr_ > 0)) FatalErrorInFunction << "Pr must be positive" << exit(FatalError);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandConstPrTransport<Thermo>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry("As", As_);
        os.writeEntry("Ts", Ts_);
        os.writeEntry("Pr", Pr_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandConstPrTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
