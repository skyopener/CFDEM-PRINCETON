Info<< "Entering UEqn" << endl;
	       
    // Gas momentum equation
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
//	/fvc::interpolate(voidfraction + scalar(0.001));

    // Source term
    fsource.internalField() = particleCloud.momCoupleM(1).expMomSource();

    //if(momentumCorrectionPeriodicFlow)
    //{		
	#include "momCorr.H"
    //}

    fsource.correctBoundaryConditions();

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
//      + (fvc::grad(voidfraction)/(fvc::average(voidfraction) + scalar(0.001)) & Rca)
      ==
      - fvm::Sp(Ksl/rho,U)
      - fsource/rho
    );

//  g                          		    // Buoyancy term transfered to p-equation    
//  Ksl/rho*Us			            // Explicit drag transfered to p-equation   

    Info << "fsource in UEqns.H = " << fsource[0] << endl;

    // Relaxing only if not performing the final iteration
    if (finalIter)
    {
        UEqn().relax(1);
    }
    else
    {
        UEqn().relax();
    }
    
    // Solving momentum predictor using flux reconstruction
    if (momentumPredictor)
    {
	solve
	(
	    UEqn() 
	 == 
	    fvc::reconstruct
	    (
	     (  g & mesh.Sf() )
	      - fvc::snGrad(p)*mesh.magSf()
	     ),
            mesh.solver(U.select(finalIter))  // Selecting linear solver settings 
	);

     }       
Info<< "Leaving UEqn" << endl;
    
