                volScalarField rUA = 1.0/UEqn().A();

                surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);

                U = rUA*UEqn().H();

                if (nCorrSoph <= 1)
                {
                    UEqn.clear();
                }

                phi = (fvc::interpolate(U*voidfraction) & mesh.Sf() )
                     + fvc::ddtPhiCorr(rUAvoidfraction, U, phi);

                surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS 
                                                - (fvc::interpolate(voidfraction*fsource/rho) & mesh.Sf()) );

                rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + particleCloud.ddtVoidfraction()
                    );
                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phiGes -= pEqn.flux();
                    }

                } // end non-orthogonal corrector loop

                #include "continuityErrorPhiPU.H"

                // Explicitly relax pressure for momentum corrector
                p.relax();

                U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA + voidfraction*fsource/rho*rUA;

                U.correctBoundaryConditions();
