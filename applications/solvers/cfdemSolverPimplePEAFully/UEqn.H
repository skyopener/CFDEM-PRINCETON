Info<< "Entering UEqn" << endl;

     fvVectorMatrix UEqn
     (
	 fvm::ddt(voidfraction, U) + fvm::div(voidfractionPhi, U)
       - fvm::Sp(fvc::ddt(voidfraction) + fvc::div(voidfractionPhi), U)
       + particleCloud.divVoidfractionTau(U, voidfraction)
//       + continuousPhaseTurbulence->divDevRhoReff(U)
      ==
       - fvm::Sp(Ksl/rho,U) 
     );

     UEqn.relax();

     volScalarField rAU(1.0/UEqn.A());
     surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

     surfaceScalarField phiForces
     (
	rAUf*(g & mesh.Sf())
      +(fvc::interpolate(rAU*Ksl/rho*Us) & mesh.Sf())
      -(fvc::interpolate(rAU*fsource/rho)& mesh.Sf())
     );

     if (pimple.momentumPredictor())
     {
	 solve
	 (
             UEqn
	  ==
             fvc::reconstruct
             (
        	 phiForces/rAUf - fvc::snGrad(p)*mesh.magSf()
             )
	 );
     }
	                
Info<< "Leaving UEqn" << endl;
    
