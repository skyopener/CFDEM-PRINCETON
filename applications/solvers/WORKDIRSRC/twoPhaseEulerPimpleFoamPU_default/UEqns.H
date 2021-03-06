    // Phase a momentum equation
    volTensorField gradUaT = fvc::grad(Ua)().T();

    //_AO_09/08/2014 - if (kineticTheory.on())
    if (kineticTheory.on() || kineticTheory.PU())
    {
	kineticTheory.solve(gradUaT);
	nuEffa = kineticTheory.mua()/rhoa;
    }
    else // If not using kinetic theory is using Ct model
    {
	nuEffa = sqr(Ct)*nutb + nua;
    }

    volTensorField Rca
    (
	"Rca",
	((2.0/3.0)*I)*(sqr(Ct)*k + nuEffa*tr(gradUaT)) - nuEffa*gradUaT
    );

    // Adding bulk viscosity if kinetic theory is used
    if (kineticTheory.on())
    {
	Rca -= ((kineticTheory.lambda()/rhoa)*tr(gradUaT))*tensor(I);
    }
    //_AO_09/03/2014 Garzo-Dufty kinetic theory OR modified kinetic theory by Chialvo-Sundaresan
    else if(kineticTheory.PU())
    {
	Rca = tensor(I) & ( kineticTheory.tau() / rhoa );
    }

    surfaceScalarField phiRa =
	-fvc::interpolate(nuEffa)*mesh.magSf()*fvc::snGrad(alpha)
	/fvc::interpolate(alpha + scalar(0.001));

    tmp<fvVectorMatrix> UaEqn
    (
	(scalar(1) + Cvm*rhob*beta/rhoa)*
	(
	    fvm::ddt(Ua)
	  + fvm::div(phia, Ua, "div(phia,Ua)")
	  - fvm::Sp(fvc::div(phia), Ua)
	)

      - fvm::laplacian(nuEffa, Ua)
      + fvc::div(Rca)

      + fvm::div(phiRa, Ua, "div(phia,Ua)")
      - fvm::Sp(fvc::div(phiRa), Ua)
      + (fvc::grad(alpha)/(fvc::average(alpha) + scalar(0.001)) & Rca)
      ==
    //  g                          // Buoyancy term transfered to p-equation
      - fvm::Sp(beta/rhoa*K, Ua)
    //+ beta/rhoa*K*Ub             // Explicit drag transfered to p-equation
      - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)
    );

    // Phase b momentum equation
    volTensorField gradUbT = fvc::grad(Ub)().T();
    volTensorField Rcb
    (
	"Rcb",
	((2.0/3.0)*I)*(k + nuEffb*tr(gradUbT)) - nuEffb*gradUbT
    );

    surfaceScalarField phiRb =
	-fvc::interpolate(nuEffb)*mesh.magSf()*fvc::snGrad(beta)
	/fvc::interpolate(beta + scalar(0.001));

    tmp<fvVectorMatrix> UbEqn
    (
	(scalar(1) + Cvm*rhob*alpha/rhob)*
	(
	    fvm::ddt(Ub)
	  + fvm::div(phib, Ub, "div(phib,Ub)")
	  - fvm::Sp(fvc::div(phib), Ub)
	)

      - fvm::laplacian(nuEffb, Ub)
      + fvc::div(Rcb)

      + fvm::div(phiRb, Ub, "div(phib,Ub)")
      - fvm::Sp(fvc::div(phiRb), Ub)

      + (fvc::grad(beta)/(fvc::average(beta) + scalar(0.001)) & Rcb)
      ==
    //  g                          // Buoyancy term transfered to p-equation
      - fvm::Sp(alpha/rhob*K, Ub)
    //+ alpha/rhob*K*Ua            // Explicit drag transfered to p-equation
      + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
    );

    // Relaxing only if not performing the final iteration
    if (finalIter)
    {
        UaEqn().relax(1);
	UbEqn().relax(1);
    }
    else
    {
        UaEqn().relax();
	UbEqn().relax();
    }
    
    // Solving momentum predictor using flux reconstruction
    if (momentumPredictor)
    {
	solve
	(
	    UaEqn() 
	 == 
	    fvc::reconstruct
	    (
		(g & mesh.Sf())
	      - fvc::snGrad(p)*mesh.magSf()/rhoa
	      - ppMagfUEqn*fvc::snGrad(alpha)*mesh.magSf()
	      + fvc::interpolate(beta/rhoa*K)*phib 
	    ),
            mesh.solver(Ua.select(finalIter))  // Selecting linear solver settings 
	);

	solve
	(
	    UbEqn() 
	 == 
	    fvc::reconstruct
	    (
		(g & mesh.Sf())
	      - fvc::snGrad(p)*mesh.magSf()/rhob
	      +	fvc::interpolate(alpha/rhob*K)*phia
	    ),
	    mesh.solver(Ub.select(finalIter)) // Selecting linear solver settings 
	);
    }
