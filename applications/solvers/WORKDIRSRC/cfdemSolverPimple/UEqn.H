Info<< "Entering UEqn" << endl;

    // Momentum predictor
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(voidfraction,U)
      + fvm::div(phi, U)
//              + turbulence->divDevReff(U)
      + particleCloud.divVoidfractionTau(U, voidfraction)
     ==
      - fvm::Sp(Ksl/rho,U)
    );

    UEqn().relax();

    // Source term
    fsource.internalField() = particleCloud.momCoupleM(1).expMomSource();

    if(momentumCorrectionPeriodicFlow)
    {
        #include "momCorr.H"
    }

    fsource.correctBoundaryConditions();

    Info << "fsource in UEqns.H at cell[0]= " << fsource[0] << endl;

    UEqn() == - voidfraction*fvc::grad(p) + Ksl/rho*Us - voidfraction*fsource/rho;

    if (momentumPredictor || pimple.momentumPredictor())
    {
        solve(UEqn());
    }
	       
     
Info<< "Leaving UEqn" << endl;
    
