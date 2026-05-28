/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "flameletThermo.H"
#include "pdfFlameletThermo.H"
#include "pureMixture.H"

#include "forGases.H"
#include "makeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // OpenFOAM 13 makeThermo uses BaseThermo::DerivedThermoType, so
    // pdfFlameletThermo is selected through flameletThermo::DerivedThermoType.
    forGases(makeThermo, flameletThermo, pureMixture);
}

// ************************************************************************* //
