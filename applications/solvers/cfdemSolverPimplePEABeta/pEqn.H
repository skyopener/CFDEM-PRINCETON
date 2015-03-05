Info<< "Entering pEqn" << endl;

    //surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

    volScalarField rUA = 1.0/UEqn().A();

    phi == (fvc::interpolate(U*voidfraction) & mesh.Sf());

    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();
    
    // Clearing temporary momentum predictor objects
    if (nCorr <= 1)
    {
	UEqn.clear();
    }

    surfaceScalarField phiDrag =
        fvc::interpolate(Ksl/rho*rUA)*(fvc::interpolate(Us*voidfraction) & mesh.Sf()) 
     - rUAf*( fvc::interpolate(fsource*voidfraction/rho) & mesh.Sf() )
     + rUAf*( g & mesh.Sf() ); 
	
    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDrag.boundaryField()[patchi] = 0.0;
        }
    }

    phi = (fvc::interpolate(U*voidfraction) & mesh.Sf()) 
         + fvc::ddtPhiCorr(rUA*voidfraction, U, phi) 
	 + phiDrag ;
    
//    surfaceScalarField Dp("(rho*(1|A(U)))", voidfractionf*voidfractionf*rUAf );
    surfaceScalarField Dp("(rho*(1|A(U)))", voidfractionf*rUAf );

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi) + particleCloud.ddtVoidfraction()
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
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phi -= rUAf*SfGradp;
            
	    p.relax();

            U += fvc::reconstruct(phiDrag - rUAf*SfGradp);
            U.correctBoundaryConditions();

        } 
    } 
	 
Info<< "Leaving pEqn" << endl;
