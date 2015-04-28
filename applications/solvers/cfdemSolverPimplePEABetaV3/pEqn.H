Info<< "Entering pEqn" << endl;

    // H/A
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    // Clearing temporary momentum predictor objects
    if (nCorr <= 1)
    {
       UEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           (fvc::interpolate(HbyA) & mesh.Sf())
         + fvc::ddtPhiCorr(rAU*voidfraction, U, phi)
         + phiForces
        )
    );
	
    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiForces.boundaryField()[patchi] = 0.0;
        }
    }

//    surfaceScalarField Dp("(rho*(1|A(U)))", voidfractionf*voidfractionf*rAUf );
    surfaceScalarField Dp("(voidfractionf*(1|A(U)))", voidfractionf*rAUf );

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(voidfractionf*phi) 
//				     + fvc::ddt(voidfraction) 	
				     + particleCloud.ddtVoidfraction()
//            fvm::laplacian(Dp, p) == fvc::div(phi) 
        );

        pEqn.setReference(pRefCell, pRefValue);
	
        pEqn.solve
	(
	    mesh.solver
	    (
		// Selecting linear solver
		p.select
		(
		     finalIter		 	// Final outer corrector
		  && corr == nCorr-1		// Final inner corrector
		  && nonOrth == nNonOrthCorr	// Final non-orhogonal corrector
		)
	    )
	);

        if (nonOrth == nNonOrthCorr)
        {
            
            phi = phiHbyA - pEqn.flux()/voidfractionf;

            p.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phiForces - pEqn.flux()/voidfractionf)/rAUf);

            U.correctBoundaryConditions();

        } 
    } 
	 
Info<< "Leaving pEqn" << endl;
