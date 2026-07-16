#include "SLFMFluid.H"

#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcVolumeIntegrate.H"
#include "surfaceInterpolate.H"

#include "fvmDiv.H"
#include "fvmLaplacian.H"

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(SLFMFluid, 0);
    addToRunTimeSelectionTable(solver, SLFMFluid, fvMesh);
}
}


Foam::solvers::SLFMFluid::SLFMFluid(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(flameletThermo::New(mesh).ptr())
    ),

    flameletThermo_(refCast<flameletThermo>(isothermalFluid::thermo_)),

    Z_(flameletThermo_.Z()),
    Zvar_(flameletThermo_.Zvar()),
    chi_st_(flameletThermo_.chi_st()),
    H_(flameletThermo_.H()),
    as_(flameletThermo_.as()),
    psi_(flameletThermo_.psi()),

    flameletsProperties_
    (
        IOobject
        (
            "flameletsProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    monitoring_(flameletsProperties_.lookup("monitoring")),
    Cg_(readScalar(flameletsProperties_.lookup("Cg"))),
    Cd_(readScalar(flameletsProperties_.lookup("Cd"))),
    Cx_(readScalar(flameletsProperties_.lookup("Cx"))),

    sigmat_
    (
        "sigmat",
        dimless,
        readScalar(flameletsProperties_.lookup("sigmat"))
    ),

    k_small_(flameletsProperties_.lookup("k_small")),
    epsilon_small_(flameletsProperties_.lookup("epsilon_small")),
    Tenv_(flameletsProperties_.lookup("Tenv")),
    varianceEquation_(flameletsProperties_.lookup("varianceEquation")),
    adiabaticMode_(flameletsProperties_.lookup("adiabaticMode")),
    radiationMode_(flameletsProperties_.lookup("radiationMode")),
    switchQoff_(false),

    Qrad_
    (
        IOobject
        (
            "Qrad",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "qradiation",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0.0
        )
    ),

    sigmaSB_
    (
        "stefan_boltzmann",
        dimensionSet(1, 0, -3, -4, 0, 0, 0),
        5.6704e-8
    ),

    kMax_
    (
        "kMax",
        k_small_.dimensions(),
        flameletsProperties_.lookupOrDefault<scalar>("kMax", GREAT)
    ),

    epsilonMax_
    (
        "epsilonMax",
        epsilon_small_.dimensions(),
        flameletsProperties_.lookupOrDefault<scalar>("epsilonMax", GREAT)
    ),

    massFlowDict_
    (
        IOobject
        (
            "massFlowProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    inletPatches_(massFlowDict_.lookup("inletPatches")),
    outletPatches_(massFlowDict_.lookup("outletPatches")),
    outInfoPtr_(new OFstream("outInfo.case"))
{
    Info<< "Reading SLFM/flamelet properties" << nl << endl;

    flameletThermo_.validate(type(), "h", "e");

    // The base class already created rho, U, phi, p, MRF, pressure reference,
    // and the compressible momentum-transport model in the OF14-native order.
    mesh.schemes().setFluxRequired(Z_.name());
    mesh.schemes().setFluxRequired(Zvar_.name());
    mesh.schemes().setFluxRequired(H_.name());

    Info<< "SLFM turbulence clamp limits: kMax=" << kMax_.value()
        << ", epsilonMax=" << epsilonMax_.value() << nl << endl;

    OFstream& outInfo = outInfoPtr_();
    outInfo << "time\t"
            << "m_in\t"
            << "v_in\t"
            << "H_in\t"
            << "m_out\t"
            << "v_out\t"
            << "H_out\t"
            << "m_net\t"
            << "v_net\t"
            << "H_net"
            << endl;
}


Foam::solvers::SLFMFluid::~SLFMFluid()
{}


void Foam::solvers::SLFMFluid::boundTurbulenceFields(const word& where)
{
    if (!mesh.foundObject<volScalarField>("k") || !mesh.foundObject<volScalarField>("epsilon"))
    {
        return;
    }

    volScalarField& k = mesh.lookupObjectRef<volScalarField>("k");
    volScalarField& epsilon = mesh.lookupObjectRef<volScalarField>("epsilon");

    const scalar kBeforeMin = gMin(k.primitiveField());
    const scalar kBeforeMax = gMax(k.primitiveField());
    const scalar eBeforeMin = gMin(epsilon.primitiveField());
    const scalar eBeforeMax = gMax(epsilon.primitiveField());

    k.max(k_small_);
    if (kMax_.value() < GREAT/10)
    {
        k.min(kMax_);
    }

    epsilon.max(epsilon_small_);
    if (epsilonMax_.value() < GREAT/10)
    {
        epsilon.min(epsilonMax_);
    }

    Info<< "SLFM clamp " << where << " k: before ["
        << kBeforeMin << ", " << kBeforeMax << "] after ["
        << gMin(k.primitiveField()) << ", " << gMax(k.primitiveField()) << "]" << nl;

    Info<< "SLFM clamp " << where << " epsilon: before ["
        << eBeforeMin << ", " << eBeforeMax << "] after ["
        << gMin(epsilon.primitiveField()) << ", " << gMax(epsilon.primitiveField()) << "]" << endl;
}


void Foam::solvers::SLFMFluid::solveMixtureFraction()
{
    boundTurbulenceFields("before-SLFM");

    const volScalarField& rho = isothermalFluid::rho;
    const surfaceScalarField& phi = isothermalFluid::phi;
    const compressible::momentumTransportModel& turb = momentumTransport();

    fvScalarMatrix ZEqn
    (
        fvm::div(phi, Z_)
      - fvm::laplacian(rho*turb.nuEff()/sigmat_, Z_)
    );

    ZEqn.relax();
    fvConstraints().constrain(ZEqn);
    ZEqn.solve("Z");
    fvConstraints().constrain(Z_);

    if (varianceEquation_ == false)
    {
        fvScalarMatrix ZvarEqn
        (
            fvm::div(phi, Zvar_)
          - fvm::laplacian(rho*turb.nut()/sigmat_, Zvar_)
          - Cg_*rho*turb.nut()*magSqr(fvc::grad(Z_))
          + Cd_*rho*turb.epsilon()/(turb.k() + k_small_)*Zvar_
        );

        ZvarEqn.relax();
        fvConstraints().constrain(ZvarEqn);
        ZvarEqn.solve("Zvar");
        fvConstraints().constrain(Zvar_);
    }
    else
    {
        Info<< "Algebraic equation for Zvar" << endl;

        scalarField& ZvarCells = Zvar_.primitiveFieldRef();

        ZvarCells =
            Cg_*rho.primitiveField()*turb.nut()().primitiveField()
           *magSqr(fvc::grad(Z_)())().primitiveField()
           *turb.k()().primitiveField()
          /(
                Cd_*rho.primitiveField()
               *(turb.epsilon()().primitiveField() + epsilon_small_.value())
           );
    }

    chi_st_ = Cx_*turb.epsilon()/(turb.k() + k_small_)*Zvar_;
}


void Foam::solvers::SLFMFluid::solveEnthalpy()
{
    const volScalarField& rho = isothermalFluid::rho;
    const surfaceScalarField& phi = isothermalFluid::phi;
    const compressible::momentumTransportModel& turb = momentumTransport();

    if (adiabaticMode_ == true)
    {
        fvScalarMatrix hEqn
        (
            fvm::div(phi, H_)
          - fvm::laplacian(rho*turb.nuEff()/sigmat_, H_)
        );

        hEqn.relax();
        fvConstraints().constrain(hEqn);
        hEqn.solve("H");
        fvConstraints().constrain(H_);
    }
    else
    {
        if (radiationMode_ == true)
        {
            Qrad_ = 4.0*sigmaSB_*as_*(pow(flameletThermo_.T(), 4) - pow(Tenv_, 4));
        }
        else if (radiationMode_ == false && switchQoff_ == false)
        {
            Qrad_ = 0.0*4.0*sigmaSB_*as_*(pow(flameletThermo_.T(), 4) - pow(Tenv_, 4));
            switchQoff_ = true;
        }

        fvScalarMatrix hEqn
        (
            fvm::div(phi, H_)
          - fvm::laplacian(rho*turb.nuEff()/sigmat_, H_)
          + Qrad_
        );

        hEqn.relax();
        fvConstraints().constrain(hEqn);
        hEqn.solve("H");
        fvConstraints().constrain(H_);
    }

    flameletThermo_.correct();
}


void Foam::solvers::SLFMFluid::thermophysicalPredictor()
{
    solveMixtureFraction();
    solveEnthalpy();
}


void Foam::solvers::SLFMFluid::momentumTransportCorrector()
{
    isothermalFluid::momentumTransportCorrector();
    boundTurbulenceFields("after-momentumTransportCorrector");
}


void Foam::solvers::SLFMFluid::writeMassFlow()
{
    scalar m_in = 0;
    scalar v_in = 0;
    scalar H_in = 0;
    scalar m_out = 0;
    scalar v_out = 0;
    scalar H_out = 0;

    const volScalarField& rho = isothermalFluid::rho;
    const volVectorField& U = isothermalFluid::U;

    OFstream& outInfo = outInfoPtr_();
    outInfo << runTime.value() << "\t";

    forAll(mesh.boundary(), patchi)
    {
        forAll(inletPatches_, iP)
        {
            if (mesh.boundary()[patchi].name() == inletPatches_[iP])
            {
                scalar flow = sum
                (
                    rho.boundaryField()[patchi]
                   *(U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf())
                );

                scalar Hflow = sum
                (
                    rho.boundaryField()[patchi]
                   *H_.boundaryField()[patchi]
                   *(U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf())
                );

                scalar volFlow = sum
                (
                    U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf()
                );

                m_in += mag(flow);
                v_in += mag(volFlow);
                H_in += mag(Hflow);
            }
        }

        forAll(outletPatches_, oP)
        {
            if (mesh.boundary()[patchi].name() == outletPatches_[oP])
            {
                scalar flow = sum
                (
                    rho.boundaryField()[patchi]
                   *(U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf())
                );

                scalar Hflow = sum
                (
                    rho.boundaryField()[patchi]
                   *H_.boundaryField()[patchi]
                   *(U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf())
                );

                scalar volFlow = sum
                (
                    U.boundaryField()[patchi] & mesh.boundary()[patchi].Sf()
                );

                m_out += mag(flow);
                v_out += mag(volFlow);
                H_out += mag(Hflow);
            }
        }
    }

    outInfo << m_in << "\t" << v_in << "\t" << H_in << "\t";
    outInfo << m_out << "\t" << v_out << "\t" << H_out << "\t";
    outInfo << m_out - m_in << "\t" << v_out - v_in << "\t" << H_out - H_in << endl;
}


void Foam::solvers::SLFMFluid::postSolve()
{
    isothermalFluid::postSolve();
    writeMassFlow();
}

// ************************************************************************* //
