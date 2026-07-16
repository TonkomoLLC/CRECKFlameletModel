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
#include "makeFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // OpenFOAM 14 fluid thermos must be registered in the basicThermo,
    // fluidThermo and flameletThermo run-time selection tables.
    // pdfFlameletThermo is selected through
    // flameletThermo::DerivedThermoType.
    forGases(makeFluidThermo, flameletThermo, pureMixture);
}

// ************************************************************************* //
