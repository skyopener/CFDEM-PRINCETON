{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic = phi;
    surfaceScalarField phir = phia - phib;

    if (g0.value() > 0.0 || kineticTheory.on())
    {
        surfaceScalarField alphaf = fvc::interpolate(alpha);
        surfaceScalarField phipp = ppMagf*fvc::snGrad(alpha)*mesh.magSf();
        phir += phipp;
        phic += fvc::interpolate(alpha)*phipp;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phic, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, beta, schemer), alpha, schemer)
        );

        if ((g0.value() > 0.0 && !kineticTheory.on()))
        {
            ppMagf = rUaAf*fvc::interpolate
            (
                (1.0/(rhoa*(alpha + scalar(0.0001))))
               *g0*min(exp(preAlphaExp*(alpha - alphaMax)), expMax)
            );
	    
	    ppMagfUEqn = fvc::interpolate
            (
                (1.0/(rhoa*(alpha + scalar(0.0001))))
               *g0*min(exp(preAlphaExp*(alpha - alphaMax)), expMax)
            );

            alphaEqn -= fvm::laplacian
            (
                (fvc::interpolate(alpha) + scalar(0.0001))*ppMagf,
                alpha,
                "laplacian(alphaPpMag,alpha)"
            );
        }
        else if (kineticTheory.on())
	{
	    ppMagf = rUaAf*fvc::interpolate
	    (
		(1.0/(rhoa*(alpha + scalar(0.0001))))
	       *kineticTheory.ppMagf(alpha)
	    );
	    
	    ppMagfUEqn = fvc::interpolate
	    (
		(1.0/(rhoa*(alpha + scalar(0.0001))))
	       *kineticTheory.ppMagf(alpha)
	    );
	    
            alphaEqn -= fvm::laplacian
            (
                (fvc::interpolate(alpha) + scalar(0.0001))*ppMagf,
                alpha,
                "laplacian(alphaPpMag,alpha)"
            );
	}

        // Relaxing only if not performing the final iteration
	if (finalIter)
	{
	    alphaEqn.relax(1);
	}
	else
	{
	    alphaEqn.relax();
	}
	
	alphaEqn.solve
	(
	    // Selecting linear solver settings  
	    mesh.solver(alpha.select(finalIter))
	);

        beta = scalar(1) - alpha;

        Info<< "Dispersed phase volume fraction = "
            << alpha.weightedAverage(mesh.V()).value()
            << "  Min(alpha) = " << min(alpha).value()
            << "  Max(alpha) = " << max(alpha).value()
            << endl;
    }
}

rho = alpha*rhoa + beta*rhob;
