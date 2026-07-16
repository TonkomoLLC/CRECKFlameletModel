/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is distributed under the GNU General Public License.

Class
    Foam::solvers::SLFMFluid

Description
    Runtime-selectable foamRun SLFM module built on the OpenFOAM-14
    isothermalFluid solver infrastructure.

\*---------------------------------------------------------------------------*/

#ifndef SLFMFluid_H
#define SLFMFluid_H

#include "isothermalFluid.H"
#include "flameletThermo.H"
#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

namespace Foam
{
namespace solvers
{

class SLFMFluid
:
    public isothermalFluid
{
protected:

    //- Flamelet thermo, down-cast from isothermalFluid::thermo_
    flameletThermo& flameletThermo_;

    //- Flamelet transport fields
    volScalarField& Z_;
    volScalarField& Zvar_;
    volScalarField& chi_st_;
    volScalarField& H_;
    volScalarField& as_;
    const volScalarField& psi_;

    //- SLFM controls
    IOdictionary flameletsProperties_;

    Switch monitoring_;
    scalar Cg_;
    scalar Cd_;
    scalar Cx_;
    dimensionedScalar sigmat_;
    dimensionedScalar k_small_;
    dimensionedScalar epsilon_small_;
    dimensionedScalar Tenv_;
    Switch varianceEquation_;
    Switch adiabaticMode_;
    Switch radiationMode_;
    bool switchQoff_;

    volScalarField Qrad_;
    dimensionedScalar sigmaSB_;

    //- Optional diagnostic turbulence clamps
    dimensionedScalar kMax_;
    dimensionedScalar epsilonMax_;

    //- Mass-flow reporting
    IOdictionary massFlowDict_;
    wordList inletPatches_;
    wordList outletPatches_;
    autoPtr<OFstream> outInfoPtr_;

    // Private Member Functions
    void boundTurbulenceFields(const word& where);
    void solveMixtureFraction();
    void solveEnthalpy();
    void writeMassFlow();

public:

    TypeName("SLFMFluid");

    SLFMFluid(fvMesh& mesh);

    SLFMFluid(const SLFMFluid&) = delete;

    virtual ~SLFMFluid();

    void operator=(const SLFMFluid&) = delete;

    // Use isothermalFluid for preSolve, momentum, pressure and turbulence.
    virtual void thermophysicalPredictor();
    virtual void momentumTransportCorrector();
    virtual void postSolve();
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
