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
    // OpenFOAM-10 makeThermo uses BaseThermo::composite as the first
    // template argument of CThermo.  flameletThermo therefore defines
    // flameletThermo::composite in flameletThermo.H.
    //
    // Keep this to gas thermo packages.  The flamelet table supplies the
    // actual thermodynamic/transport fields; liquid/tabulated/polynomial
    // OpenFOAM thermo packages are unnecessary here and greatly expand the
    // template instantiation set.
    forGases(makeThermo, flameletThermo, pdfFlameletThermo, pureMixture);
}

// ************************************************************************* //
