/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postCalc

Description
    Generic wrapper for calculating a quantity at each time

\*---------------------------------------------------------------------------*/

#include "calcCollisionalForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	

	Foam::timeSelector::addOptions();
	#   include "addRegionOption.H"
	Foam::argList::addBoolOption
	(
	"noWrite",
	"suppress writing results"
	);
	#include "addDictOption.H"

	#include "setRootCase.H"
	#include "createTime.H"
	Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

	#include "createFields.H" 	 

	// Create particle cloud			
	cfdemCloud particleCloud(mesh);
	
	particleCloud.reAllocArrays();

	double **positions;
	double **velocities;
	double **omegas;
	double **radii;
	double **cellID;
	
	double **forces;
	double **types;

	particleCloud.dataExchangeM().allocateArray(positions,0.,3);
	particleCloud.dataExchangeM().allocateArray(velocities,0.,3);
	particleCloud.dataExchangeM().allocateArray(omegas,0.,3);
	particleCloud.get_radii(radii); 
	particleCloud.get_cellIDs(cellID);
	
	particleCloud.dataExchangeM().allocateArray(forces,0.,3);
	particleCloud.dataExchangeM().allocateArray(types,0.,1);
				
	// Read dictionary + dictionaryProps 
	const dictionary dict(particleCloud.couplingProperties());
	int maxNumberOfParticles(0);
	
	// oneWayVTK dictionary
	const dictionary oneWayVTKPropsDict(dict.subDict("oneWayVTKProps"));
	maxNumberOfParticles=readScalar(oneWayVTKPropsDict.lookup("maxNumberOfParticles"));
	
	// Read particleStress Sub-dictionary
	const dictionary ParticleStressPropsDict(dict.subDict("ParticleStressProps"));
	const scalar rhop(readScalar(ParticleStressPropsDict.lookup("rhoParticle")));					
	
	// Debuging
	bool verbose(false);
	// Particle ID for debuging
	int exIndex(0);
	if(ParticleStressPropsDict.found("verbose")) 
	{
		verbose = true;	
		exIndex = readInt(ParticleStressPropsDict.lookup("exIndex"));	
	}	
	
	// Collision part debugging
	bool verboseColl(false);
	if(ParticleStressPropsDict.found("verboseColl")) 
	{
		verboseColl = true;			
	}	
			
	// Number particle class
	int nParticleClass = 1;
	bool bidisperse(false);
	// Bidisperse case
	if(ParticleStressPropsDict.found("nParticleClass"))
	{
		nParticleClass = readInt(ParticleStressPropsDict.lookup("nParticleClass"));
		bidisperse = true;
		Info << " " << endl;
		Info << "Bi-disperse case, number of particle classes = " << nParticleClass << endl;		
	}
	
	// Collision dictionary
	scalar collisionDp(1);
	bool calcCollision(false);
	if(ParticleStressPropsDict.found("calcCollision"))
	{
		calcCollision=true;
		collisionDp=readScalar(ParticleStressPropsDict.lookup("collisionDp"));
	}

	//
	const dictionary twoWayMPIPropsDict(dict.subDict("twoWayMPIProps"));	
	// DEM input file parameters
	fileName DEMinputFilename(twoWayMPIPropsDict.lookup("liggghtsPath"));

	// Create output folder	
        fileName outputRelativePath(ParticleStressPropsDict.lookup("outputRelativePath"));
	if( !isDir(mesh.time().path()/outputRelativePath) )
	{
		mkDir(mesh.time().path()/outputRelativePath );
	}

    	OFstream* outputKinFile;
	OFstream* outputCollFile;

    	outputKinFile =  new OFstream(mesh.time().path()/outputRelativePath/"kineticStress");
	*outputKinFile  << "#Time  " << tab << " " 
	                << "Kin_XX " << tab << " "
			<< "Kin_YY " << tab << " "
			<< "Kin_ZZ " << tab << " " 
			<< "Kin_XY " << tab << " " 
			<< "Kin_XZ " << tab << " " 
			<< "Kin_YZ " << tab << " "; 		

	if( bidisperse ) 
	{
		for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        {
			*outputKinFile  << "Kin_XX["<<iPartClass<<"] " << tab << " "
					<< "Kin_YY["<<iPartClass<<"] " << tab << " "
					<< "Kin_ZZ["<<iPartClass<<"] " << tab << " " 
					<< "Kin_XY["<<iPartClass<<"] " << tab << " " 
					<< "Kin_XZ["<<iPartClass<<"] " << tab << " " 
					<< "Kin_YZ["<<iPartClass<<"] " << tab << " "; 
		}
	}
	*outputKinFile << endl;
	
	if( calcCollision )
	{
    		outputCollFile =  new OFstream(mesh.time().path()/outputRelativePath/"collisonalStress");
		*outputCollFile  << "#Time  " << tab << " " 
	                	 << "Coll_XX " << tab << " "
				 << "Coll_YY " << tab << " "
				 << "Coll_ZZ " << tab << " " 
				 << "Coll_XY " << tab << " " 
				 << "Coll_XZ " << tab << " " 
				 << "Coll_YZ " << tab << " "; 	
		if( bidisperse ) 
		{
			for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        	{
				*outputCollFile  << "Coll_XX["<<iPartClass<<"] " << tab << " "
						 << "Coll_YY["<<iPartClass<<"] " << tab << " "
						 << "Coll_ZZ["<<iPartClass<<"] " << tab << " " 
						 << "Coll_XY["<<iPartClass<<"] " << tab << " " 
						 << "Coll_XZ["<<iPartClass<<"] " << tab << " " 
						 << "Coll_YZ["<<iPartClass<<"] " << tab << " "; 			
			}
		}
		*outputCollFile << endl;	
	}
		
	// Collision parameters
	
	double* youngsModulus = new double[nParticleClass];
	double* poissonsRatio = new double[nParticleClass];
		
	double **coefficientRestitution = new double*[nParticleClass];
	for(int i = 0; i < nParticleClass; ++i) 
	{
    		coefficientRestitution[i] = new double[nParticleClass];
	}	

	double **coefficientFriction = new double*[nParticleClass];
	for(int i = 0; i < nParticleClass; ++i) 
	{
    		coefficientFriction[i] = new double[nParticleClass];
	}
		
	//double* coefficientRestitution = new double[nParticleClass*nParticleClass]; 
	//double* coefficientFriction = new double[nParticleClass];
	
	double k_n(0); // = new double[nParticleClass];
	double k_t(0); // = new double[nParticleClass];
	double gamma_n(0); // = new double[nParticleClass];
	double gamma_t(0); // = new double[nParticleClass];
	double mu_f(0); //= new double[nParticleClass];
	double e_n(0); // = new double[nParticleClass];
	double e_t(0); // = new double[nParticleClass];

	double dt(0);

	bool tangential_history=false;
	if( calcCollision )
	{
		if(ParticleStressPropsDict.found("tangential_history")) 	
		{
			tangential_history = true;
			Info << "Tangential history effect is active " << endl;
		}
	}
	bool liquid_transfer=false;
	if(ParticleStressPropsDict.found("liquid_transfer")) liquid_transfer = true;	
	double surf_tension(0);
	double fluid_visc(0);

	bool cohesion=false;
	if(ParticleStressPropsDict.found("cohesion")) cohesion = true;	
	double minimumDistanceVdW(0);
	//double* cohEnergyDens = new double[nParticleClass];
	
	double **cohEnergyDens = new double*[nParticleClass];
	for(int i = 0; i < nParticleClass; ++i) 
	{
    		cohEnergyDens[i] = new double[nParticleClass];
	}
				
	// Define & Init collision variables
	#include "init_collision_variables.H" 
				
	// Collision model
	int collisionModelI(0); // Hertz (default)	
	
	if(calcCollision)
	{
		#include "collision_parameters.H"
	}
		
	// Neighboring list parameters				
	int k(0); 

	// Dimensions, exact OR approximate  
	int dim(0); double eps(0);

	// Number of points
	int nPts(0);

	// Data points
	ANNpointArray dataPts;
	// Query points
	ANNpoint queryPt;

	ANNidxArray	nnIdx;			// 	near neighbour indices
	ANNdistArray 	dists;			//	near neighbour distances
	ANNkd_tree* 	kdTree;			//	search structure
	
	// Total volume of domain
	scalar domainVol(0);
	forAll(mesh.C(),cellI)
	{
		domainVol +=mesh.V()[cellI];
	}
	Info << " " << endl;
	Info << "Domain volume[m^3] = " << domainVol << endl;	
	
	// Results variables
	// Collisional stress tensor
	symmTensor sigma_coll_JI(0,0,0,0,0,0);
	symmTensor sigma_coll[nParticleClass+1];
	// Kinetic stress tensor
	symmTensor sigma_kin[nParticleClass+1]; // nParticleClass+1 --> for mixture velocity
	
	// Mean velocities
	vector meanVel[nParticleClass+1];	// nParticleClass+1 --> for mixture velocity
	int npPartClass[nParticleClass+1];
									
	forAll(timeDirs, timeI)
	{
  
		runTime.setTime(timeDirs[timeI], timeI);

		Foam::Info << " " << endl;
		Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

		mesh.readUpdate();

		if ( runTime.timeName() != "0" )
		{
			
			int count = runTime.value() / particleCloud.dataExchangeM().DEMts();
			dt = particleCloud.dataExchangeM().DEMts(); 				

			Info<< " " << endl;
			particleCloud.dataExchangeM().getData("v","vector-atom",velocities,count);
			Info<< " 	Reading particle velocities" << endl;
			Info<< " " << endl;

			particleCloud.dataExchangeM().getData("omega","vector-atom",omegas,count);
			Info<< " 	Reading particle angular velocities" << endl;
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("x","vector-atom",positions,count);
			Info<< " 	Reading particle positions" << endl;
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("radius","scalar-atom",radii,count);
			Info<< " 	Reading particle radius" << endl;		
			Info<< " " << endl;
			
			particleCloud.dataExchangeM().getData("type","scalar-atom",types,count);
			Info<< " 	Reading particle types " << endl;
			Info<< " " << endl;

			if ( verbose )
			{
				particleCloud.dataExchangeM().getData("f","vector-atom",forces,count);
				Info<< " 	Reading forces on particles" << endl;
				Info<< " " << endl;
			}	
					
			particleCloud.locateM().findCell(NULL,positions,cellID,particleCloud.numberOfParticles());
			particleCloud.setPos(positions);
			particleCloud.setVel(velocities);
			particleCloud.setOmega(omegas);
			particleCloud.setType(types);
			
			if(verbose)
			{
				Info << "" << endl;
				Info << " index  = " << exIndex << endl;
				Info << " rp     = " << particleCloud.radius(exIndex) << endl;
				Info << " Type   = " << particleCloud.type(exIndex) << endl;
				Info << " Vp     = " << particleCloud.velocity(exIndex) << endl;
				Info << " Omegap = " << particleCloud.omega(exIndex) << endl;				
				Info << " Xp     = " << particleCloud.position(exIndex) << endl;
				Info << " CellID = " << particleCloud.particleCell(exIndex) << endl;
				Info << " Fp     = " << forces[exIndex][0] << " " << forces[exIndex][1] << " " << forces[exIndex][2] << endl;
			}				
			
			if(calcCollision)
			{
				// Neighboring list parameters
				k = 20; 
				if(particleCloud.numberOfParticles()<k) k = particleCloud.numberOfParticles();

				// Dimensions, exact OR approximate  
				dim=3; eps = 0;

				// Number of points
				nPts =  particleCloud.numberOfParticles();		       		

				// Allocate 
				queryPt = annAllocPt(dim);
				dataPts = annAllocPts(nPts, dim);
				nnIdx = new ANNidx[k];
				dists = new ANNdist[k];

				// Particle collisional stresses
				sigma_coll_JI = symmTensor(0,0,0,0,0,0);				
				for(int iPartClass = 0; iPartClass <= nParticleClass; iPartClass++)
	        		{
					sigma_coll[iPartClass] = symmTensor(0,0,0,0,0,0);
				}	
			}			
						
			// Initiate
			int iPartClass = 0;
			npPartClass[iPartClass] = particleCloud.numberOfParticles();
			sigma_kin[iPartClass] = symmTensor(0,0,0,0,0,0); 
			
			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        		{
					// Initiate particle number of each class
					npPartClass[iPartClass] = 0;
					
					// Initiate velocities
					for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] = 0 ;

					// Initiate particle kinetic stresses
					sigma_kin[iPartClass] = symmTensor(0,0,0,0,0,0); 
				}
			}		
				
			// Create particle list
			for(int index = 0; index <  particleCloud.numberOfParticles(); index++)
	        	{							

				if(calcCollision)	
				{
					dataPts[index][0] = particleCloud.position(index).x();		
					dataPts[index][1] = particleCloud.position(index).y();		
					dataPts[index][2] = particleCloud.position(index).z();

					for (int dir=0;dir<3;dir++)
					{	
						fcoll[index][dir] = 0;
					 	ftan[index][dir] = 0;		
 					 	fcap[index][dir] = 0;
 						fvisc[index][dir] = 0;
					}
				}

				// Total velocity of particles
				int iPartClass = 0;
				for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
				
				if( bidisperse ) 
				{
					iPartClass = particleCloud.type(index);
						
					for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir] += particleCloud.velocity(index)[idir];
					npPartClass[iPartClass]++;	
		
				}	

			}	
	
			// Normalize
			for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        		{
					for(int idir=0; idir<3; idir++) meanVel[iPartClass][idir]/= npPartClass[iPartClass];
				}
			}				
			if(verbose)
			{
				iPartClass = 0;									
				Info << " " << endl;
				Info << " Particle class = " << iPartClass << endl;
				Info << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
				Info << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
				Info << " <u_p,z> = " << meanVel[iPartClass][2] << endl;

				if( bidisperse ) 
				{
					for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        			{
						Info << " " << endl;
						Info << " Particle class = " << iPartClass << endl;
						Info << " <u_p,x> = " << meanVel[iPartClass][0] << endl;
						Info << " <u_p,y> = " << meanVel[iPartClass][1] << endl;
						Info << " <u_p,z> = " << meanVel[iPartClass][2] << endl;								
					}
				}
								
			}

			if(calcCollision)
			{
				// Create Tree structure
				kdTree = new ANNkd_tree(dataPts, nPts, dim);
			}
			
			vector pos_i(0,0,0);
			vector pos_j(0,0,0);
			scalar dist(0);
			
			// Particle mass
			scalar mass_p(0);
			int typeJ(0);
			int typeI(0);
			label cellI;
					
			for(int index_j = 0; index_j <  particleCloud.numberOfParticles(); index_j++)		
			{
				// Cell ID
				cellI = particleCloud.cellIDs()[index_j][0];;
				
				// Particle type
				typeJ = particleCloud.type(index_j) - 1; // Just nclass starts from "0"
				
				if(cellI > -1)
				{
					if(calcCollision)
					{
						scalar sqRad = collisionDp*particleCloud.radius(index_j);

						queryPt[0] = particleCloud.position(index_j).x();
						queryPt[1] = particleCloud.position(index_j).y();
						queryPt[2] = particleCloud.position(index_j).z();

						kdTree->annkFRSearch(
			                				queryPt,			// query point					
									sqRad,				// squared radius
									k,				// number of the near neighbours to return
									nnIdx,				// nearest neighbor array
									dists,				// dist to near neighbours
									eps			);					

						//for (int index_i = 1; index_i < k; index_i++) // k=0 particle itself
						for (int index_i = 0; index_i < k; index_i++) // k=0 particle itself ???
						{

							if ( nnIdx[index_i] != index_j )
							{
						         typeI = particleCloud.type(nnIdx[index_i])-1; // Just nclass stats from "0"

									// Calculate collision forces
									calcForce(			  
											particleCloud,
										      collisionModelI, 
											      index_j,
										       nnIdx[index_i],     
												fcoll, 
												 ftan,
												  k_n, 
												  k_t, 
											      gamma_n, 
											      gamma_t, 
                                                        				youngsModulus,
                                                        				poissonsRatio,
                                        				       coefficientRestitution,
                                                				  coefficientFriction,				
											      delta_t, 
												 mu_f, 
												   dt, 
					        				   tangential_history, 
							  					  liq, 
										      liquid_transfer, 
						        				    liquidVol, 
						        				 surf_tension, 
						        				   fluid_visc, 
							        				 fcap, 
												fvisc,
						        				  first_touch,
											     cohesion,
										   minimumDistanceVdW,
						        				cohEnergyDens,
												 fcoh,
												 rhop,
												  e_n,
												  e_t,
											sigma_coll_JI,
											   bidisperse,
								        			typeJ,
												typeI,
											  verboseColl   		);


									iPartClass = 0;
                                                			sigma_coll[iPartClass][0] += sigma_coll_JI.xx();
                                                			sigma_coll[iPartClass][1] += sigma_coll_JI.yy();
                                                			sigma_coll[iPartClass][2] += sigma_coll_JI.zz();
                                                			sigma_coll[iPartClass][3] += sigma_coll_JI.xy();
                                                			sigma_coll[iPartClass][4] += sigma_coll_JI.xz();
                                                			sigma_coll[iPartClass][5] += sigma_coll_JI.yz();


									if( bidisperse )
									{							  						
                                                        			if( typeJ == typeI )
                                                        			{
                                                 
											iPartClass = typeJ + 1 ;
											sigma_coll[iPartClass][0] += sigma_coll_JI.xx();
                                                                			sigma_coll[iPartClass][1] += sigma_coll_JI.yy();
                                                                			sigma_coll[iPartClass][2] += sigma_coll_JI.zz();
                                                                			sigma_coll[iPartClass][3] += sigma_coll_JI.xy();
                                                                			sigma_coll[iPartClass][4] += sigma_coll_JI.xz();
                                                                			sigma_coll[iPartClass][5] += sigma_coll_JI.yz();
                                                        			}

									}
									
							}	

						}

					}	

					// Particle mass
					mass_p = 4./3.*rhop*constant::mathematical::pi*particleCloud.radius(index_j)*particleCloud.radius(index_j)*particleCloud.radius(index_j);	

					// Particle kinetic stress
					int iPartClass = 0;
					sigma_kin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ); 

					sigma_kin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
									   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ); 

					sigma_kin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ) 
									   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ); 

					sigma_kin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

					sigma_kin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                        	   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

					sigma_kin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
									   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );
											   
					if( bidisperse ) 
					{				
						iPartClass = typeJ + 1;

						sigma_kin[iPartClass][0] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	   * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ); 

						sigma_kin[iPartClass][1] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
										   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ); 

						sigma_kin[iPartClass][2] += mass_p * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ) 
										   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] ); 

						sigma_kin[iPartClass][3] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	   * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] );

						sigma_kin[iPartClass][4] += mass_p * ( particleCloud.velocity(index_j).x() - meanVel[iPartClass][0] ) 
					                                	   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );

						sigma_kin[iPartClass][5] += mass_p * ( particleCloud.velocity(index_j).y() - meanVel[iPartClass][1] ) 
										   * ( particleCloud.velocity(index_j).z() - meanVel[iPartClass][2] );	
										   
					}	
				}
			}
			
			// Write output'
			Info << " " << endl;
			Info << " Writing particle kinetic stresses into the file " << mesh.time().path()/outputRelativePath/"kineticStress" << endl;
			iPartClass = 0;
			*outputKinFile  << runTime.value() 		      << tab << " " 	
					<< sigma_kin[iPartClass][0]/domainVol << tab << " " 		
					<< sigma_kin[iPartClass][1]/domainVol << tab << " " 
					<< sigma_kin[iPartClass][2]/domainVol << tab << " " 		
					<< sigma_kin[iPartClass][3]/domainVol << tab << " " 
					<< sigma_kin[iPartClass][4]/domainVol << tab << " " 		
					<< sigma_kin[iPartClass][5]/domainVol << tab << " " ;								
			
			if( bidisperse ) 
			{
				for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        		{
					*outputKinFile  << sigma_kin[iPartClass][0]/domainVol << tab << " " 		
							<< sigma_kin[iPartClass][1]/domainVol << tab << " " 
							<< sigma_kin[iPartClass][2]/domainVol << tab << " " 		
							<< sigma_kin[iPartClass][3]/domainVol << tab << " " 
							<< sigma_kin[iPartClass][4]/domainVol << tab << " " 		
							<< sigma_kin[iPartClass][5]/domainVol << tab << " " ; 
				}
			}
			*outputKinFile << endl;

			if( calcCollision )
			{
    				Info<< " Writing particle kinetic stresses into the file " << mesh.time().path()/outputRelativePath/"collisionalStress" << endl;
				iPartClass = 0;
				*outputCollFile  << runTime.value() 		        << tab << " " 	
						 << sigma_coll[iPartClass][0]/domainVol << tab << " " 		
					  	 << sigma_coll[iPartClass][1]/domainVol << tab << " " 
						 << sigma_coll[iPartClass][2]/domainVol << tab << " " 		
						 << sigma_coll[iPartClass][3]/domainVol << tab << " " 
						 << sigma_coll[iPartClass][4]/domainVol << tab << " " 		
						 << sigma_coll[iPartClass][5]/domainVol << tab << " " ;
				if( bidisperse ) 
				{
					for(int iPartClass = 1; iPartClass <= nParticleClass; iPartClass++)
	        			{
						*outputCollFile  << sigma_coll[iPartClass][0]/domainVol << tab << " " 		
					  			 << sigma_coll[iPartClass][1]/domainVol << tab << " " 
								 << sigma_coll[iPartClass][2]/domainVol << tab << " " 		
								 << sigma_coll[iPartClass][3]/domainVol << tab << " " 
								 << sigma_coll[iPartClass][4]/domainVol << tab << " " 		
								 << sigma_coll[iPartClass][5]/domainVol << tab << " " ;			
					}
				}
				*outputCollFile << endl;	
			}
						
			if(verbose)
			{
				if(calcCollision)
				{
					Info << "" << endl;
				 	Info << " Fcoll = " << fcoll[exIndex][0] << " " << fcoll[exIndex][1] << " " << fcoll[exIndex][2] << endl;
				}
				if(cohesion)
				{
					Info << " Fcoh = " <<  fcoh[exIndex][0] << " " <<  fcoh[exIndex][1] << " " <<  fcoh[exIndex][2] << endl;
					Info << "" << endl;
				}	
			
				iPartClass = 0;
				Info << " " << endl;
				Info << " Particle class = " << iPartClass << " kinetic stresses" << endl;
				Info << " sigma_kin_xx= " << sigma_kin[iPartClass][0]/domainVol << endl;
				Info << " sigma_kin_yy= " << sigma_kin[iPartClass][1]/domainVol << endl;
				Info << " sigma_kin_zz= " << sigma_kin[iPartClass][2]/domainVol << endl;			
				Info << " sigma_kin_xy= " << sigma_kin[iPartClass][3]/domainVol << endl;
				Info << " sigma_kin_xz= " << sigma_kin[iPartClass][4]/domainVol << endl;
				Info << " sigma_kin_yz= " << sigma_kin[iPartClass][5]/domainVol << endl;
				
				if(calcCollision)
				{
					Info << " " << endl;
					Info << " Particle class = " << iPartClass << " collisional stresses" << endl;
					Info << " sigma_coll_xx= " << sigma_coll[iPartClass][0]/domainVol << endl;
					Info << " sigma_coll_yy= " << sigma_coll[iPartClass][1]/domainVol << endl;
					Info << " sigma_coll_zz= " << sigma_coll[iPartClass][2]/domainVol << endl;						
					Info << " sigma_coll_xy= " << sigma_coll[iPartClass][3]/domainVol << endl;
					Info << " sigma_coll_xz= " << sigma_coll[iPartClass][4]/domainVol << endl;			
					Info << " sigma_coll_yz= " << sigma_coll[iPartClass][5]/domainVol << endl;
				}						
								
	
				if( bidisperse ) 
				{
					for(int iPartClass = 1; iPartClass <= nParticleClass ; iPartClass++)
	        			{
						Info << " " << endl;
						Info << " Particle class = " << iPartClass << " kinetic stresses" << endl;
						Info << " sigma_kin_xx= " << sigma_kin[iPartClass][0]/domainVol << endl;
						Info << " sigma_kin_yy= " << sigma_kin[iPartClass][1]/domainVol << endl;
						Info << " sigma_kin_zz= " << sigma_kin[iPartClass][2]/domainVol << endl;			
						Info << " sigma_kin_xy= " << sigma_kin[iPartClass][3]/domainVol << endl;
						Info << " sigma_kin_xz= " << sigma_kin[iPartClass][4]/domainVol << endl;
						Info << " sigma_kin_yz= " << sigma_kin[iPartClass][5]/domainVol << endl;

						if(calcCollision)
						{
							Info << " " << endl;
							Info << " Particle class = " << iPartClass << " collisional stresses" << endl;
							Info << " sigma_coll_xx= " << sigma_coll[iPartClass][0]/domainVol << endl;
							Info << " sigma_coll_yy= " << sigma_coll[iPartClass][1]/domainVol << endl;
							Info << " sigma_coll_zz= " << sigma_coll[iPartClass][2]/domainVol << endl;						
							Info << " sigma_coll_xy= " << sigma_coll[iPartClass][3]/domainVol << endl;
							Info << " sigma_coll_xz= " << sigma_coll[iPartClass][4]/domainVol << endl;			
							Info << " sigma_coll_yz= " << sigma_coll[iPartClass][5]/domainVol << endl;
						}						
					}
				}


			
			}
		}	


	}
		
	particleCloud.dataExchangeM().destroy(positions,3);
	particleCloud.dataExchangeM().destroy(velocities,3);

	Foam::Info	<<  " " << Foam::endl;    
	Foam::Info	<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
        		<< nl << Foam::endl;

	Foam::Info<< "End\n" << Foam::endl;
	
	

}


// ************************************************************************* //
