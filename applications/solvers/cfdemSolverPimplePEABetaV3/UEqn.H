Info<< "Entering UEqn" << endl;
	       
    // Gas momentum equation
    surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

    volTensorField gradUT = fvc::grad(U)().T();

    volTensorField Rca
    (
	"Rca",
	((2.0/3.0)*I)*(nufField*tr(gradUT)) - nufField*gradUT
    );

    // Adding bulk viscosity if compressibility is important
    //  lambda: bulk viscosity
    //	Rca -= (lambda/rho)*tr(gradUT))*tensor(I);

    surfaceScalarField phiR =
	-fvc::interpolate(nufField)*mesh.magSf()*fvc::snGrad(voidfraction)
	/fvc::interpolate(voidfraction);

    tmp<fvVectorMatrix> UEqn
    (
	(
	    fvm::ddt(U)
	  + fvm::div(phi,U)
	  - fvm::Sp(fvc::div(phi), U)
	)

      - fvm::laplacian(nufField, U)
      + fvc::div(Rca)

      + fvm::div(phiR, U)
      - fvm::Sp(fvc::div(phiR), U)
      + (fvc::grad(voidfraction)/(fvc::average(voidfraction)) & Rca)
      ==
      - fvm::Sp(Ksl/rho/(fvc::average(voidfraction)),U)
    );

//  g                          		// Buoyancy term transfered to momentum prediction    
//  Ksl/rho*Us			        // Explicit drag transfered to momentum prediction
//  fsource				// Source transferres to momentum prediction

    // Relaxing only if not performing the final iteration
    if (finalIter)
    {
        UEqn().relax(1);
    }
    else
    {
        UEqn().relax();
    }

    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    surfaceScalarField phiForces =
//       rAUf * (  ( fvc::interpolate(Ksl/rho*Us/voidfraction) & mesh.Sf() )
//               - ( fvc::interpolate(fsource/rho/voidfraction)& mesh.Sf() )
//               + ( g & mesh.Sf() ) );
              (  ( fvc::interpolate(rAU*Ksl/rho*Us/voidfraction) & mesh.Sf() )
               - ( fvc::interpolate(rAU*fsource/rho/voidfraction)& mesh.Sf() )
               + rAUf * ( g & mesh.Sf() ) );



    // Solving momentum predictor using flux reconstruction
    if (momentumPredictor)
    {
	solve
	(
	    UEqn() 
	 == 
	    fvc::reconstruct
	    (
	      phiForces/rAUf - ( (fvc::snGrad(p)/voidfractionf) * mesh.magSf() )	
	    ),
            mesh.solver(U.select(finalIter))  // Selecting linear solver settings 
	);

    } 
        
/*     // Solving momentum predictor using flux reconstruction
    if (momentumPredictor)
    {
	solve
	(
	    UEqn() 
	 == 
	    fvc::reconstruct
	    (
	     (  g & mesh.Sf() )
	      - ( fvc::snGrad(p)/voidfractionf ) * mesh.magSf()
	      + ( ( fvc::interpolate(Ksl/rho*Us)  / voidfractionf ) & mesh.Sf() )
	      - ( ( fvc::interpolate(fsource/rho) / voidfractionf ) & mesh.Sf() )
	     ),
            mesh.solver(U.select(finalIter))  // Selecting linear solver settings 
	);

     }  */
	                
Info<< "Leaving UEqn" << endl;
    
