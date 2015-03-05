/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "cfdemCloud.H"
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"

#include "forceModel.H"

#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    
    cfdemCloud particleCloud(mesh);
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Gas viscosity
    #ifdef comp
        const volScalarField nufField = particleCloud.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud.turbulence().nu();
    #endif

    while (runTime.loop())
    {
        Info<< "\nStarting time loop\n" << endl;
            particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPIMPLEControls.H"
        #include "CourantNo.H"

        // Particle equations
        particleCloud.clockM().start(2,"Coupling");
        particleCloud.evolve(voidfraction,Us,U);

        Info << "update Ksl.internalField()" << endl;
        Ksl.oldTime().internalField() = particleCloud.momCoupleM(0).impMomSource();
        particleCloud.smoothingM().smoothen(Ksl);
        Ksl.correctBoundaryConditions();

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

	// Source term
	fsource.internalField() = particleCloud.momCoupleM(1).expMomSource();
	fsource.correctBoundaryConditions();

	// Correct momentum deficit
	if(correctCouplingError) fsource += fresidual;	
	 
	for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            bool finalIter = oCorr == nOuterCorr-1;
	    
	    if (nOuterCorr != 1)
	    {
		p.storePrevIter();
	    }

	    #include "UEqn.H"
	    
	    // --- PISO loop
	    int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);	    
            for (int corr=0; corr<nCorrSoph; corr++)	    
            //for (int corr=0; corr<nCorr; corr++)
	    {

	    Info<< "Pressure solver loop:" << corr << endl;
	    
		#include "pEqn.H"
		#include "continuityErrs.H"
		
	    }

        }

	// Residual source term
	forAll(mesh.cells(),cellI)
	{
    	    fresidual[cellI] =  Ksl[cellI] * ( Us[cellI] - U[cellI] ) 
			      + particleCloud.forceM(0).impParticleForces()[cellI] / mesh.V()[cellI];
	} 	

	fresidual.correctBoundaryConditions();
	 
	Info << "fresidual in UEqns.H = " << 
	 "[" << min(fresidual).value() << ":" 
             << max(fresidual).value() << "]" << endl;  
   
        if(fullyPeriodicFLow)
	{
		#include "momCorr.H"
	} 
    	if(mixtureVelocityCorrectionPeriodicFlow)
    	{	
		#include "mixVelCorr.H"
        } 
	       
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;
    
    return 0;
}


// ************************************************************************* //
