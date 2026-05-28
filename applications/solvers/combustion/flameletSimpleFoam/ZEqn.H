//-  Solving equation of variance of mixture fraction
fvScalarMatrix ZEqn
(
    fvm::div(phi, Z)
  - fvm::laplacian(rho*turbulence->nuEff()/sigmat, Z)
);

ZEqn.relax();

ZEqn.solve("Z");

if (varianceEquation == false)
{
    //- Solving equation of variance of mixture fraction
    fvScalarMatrix ZvarEqn
    (
    	(
    	    fvm::div(phi, Zvar)
    	  - fvm::laplacian(rho*turbulence->nut()/sigmat, Zvar)
    	  - Cg*rho*turbulence->nut()*magSqr(fvc::grad(Z))
    	  + Cd*rho*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar
    	)
    );

    ZvarEqn.relax();

    ZvarEqn.solve("Zvar");
}
else
{
    Info << "Algebraic equation for Zvar" << endl;

    scalarField& ZvarCells = Zvar.ref();

    ZvarCells =
        Cg*rho*turbulence->nut()*magSqr(fvc::grad(Z))*turbulence->k()
      / (Cd*rho*(turbulence->epsilon()+epsilon_small));
}

//- Updating stoichiometric dissipation rate
chi_st=Cx*turbulence->epsilon()/(turbulence->k()+k_small)*Zvar;
