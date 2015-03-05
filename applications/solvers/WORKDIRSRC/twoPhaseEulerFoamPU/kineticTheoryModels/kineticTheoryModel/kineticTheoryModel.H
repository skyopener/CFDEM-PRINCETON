/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phase1_;
        const volVectorField& U1_;
        const volVectorField& U2_;
        const volScalarField& alpha1_;
        const surfaceScalarField& phi1_;

        const dragModel& drag1_;

        const dimensionedScalar& rho1_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nu1_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        //- use generation == dissipation
        Switch equilibrium_;

        //_AO_09/30/2014 - use Garzo-Dufty kinetic theory
        Switch kineticTheoryGD_;

        //_AO_09/30/2014 - use modified kinetic theory by Chialvo-Sundaresan
        Switch mofidiedKineticTheoryPU_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<kineticTheoryModels::radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        autoPtr<frictionalStressModel> frictionalStressModel_;

        //- coefficient of restitution
        const dimensionedScalar e_;

        //_AO_09/02/2014 - effective restitution coefficient
        mutable dimensionedScalar eEff_;
	
        //_AO_09/02/2014 - friction coefficient
        mutable scalar muFric_;
	
        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- min value for which the frictional stresses are zero
        const dimensionedScalar alphaMinFriction_;

        //_AO_09/02/2014 - max vol.frac for inertial regime
        //mutable dimensionedScalar alphaf_;
        mutable scalar alphaf_;

        //_AO_09/02/2014 - critical volume fraction
        mutable scalar alphac_;

        //_AO_09/02/2014 - yield stress ratio
        mutable scalar upsilons_;

        //- material constant for frictional normal stress
        const dimensionedScalar Fr_;

        //- material constant for frictional normal stress
        const dimensionedScalar eta_;

        //- material constant for frictional normal stress
        const dimensionedScalar p_;

        //- angle of internal friction
        const dimensionedScalar phi_;

        //- The granular energy/temperature
        volScalarField Theta_;

        //- The granular viscosity
        volScalarField mu1_;

        //- The granular bulk viscosity
        volScalarField lambda_;

        //- The granular pressure
        volScalarField pa_;

        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;

        //_AO_09/02/2014 - The shear stress ratio
        volScalarField upsilon_;
	
	//_AO_09/02/2014 - The shear stress 
	volSymmTensorField tau_;

	//_AO_09/02/2014 - The derivative of particle pressure
	volScalarField ppMagf_;

	//- AO_09/02/2014 -The derivative of the radial distribution function
        volScalarField gs0Prime_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phase1,
            const volVectorField& U2,
            const volScalarField& alpha1,
            const dragModel& drag1
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        void solve(const volTensorField& gradU1t);

	//_AO_09/01/2014 Garzo-Dufty Model
	//void solveGD(const volTensorField& gradU1t);
		
	//_AO_09/02/2014_PU development of kinetic theory 
	//void solvePU(const volTensorField& gradU1t);
	
        bool on() const
        {
            return kineticTheory_;
        }

	//_AO_09/01/2014 Garzo-Dufty Model
        bool GD() const
        {
            return kineticTheoryGD_;
        }
	
 	//_AO_09/30/2014 modified kinetic theory		
        bool PU() const
        {
            return mofidiedKineticTheoryPU_;
        }

        const volScalarField& mu1() const
        {
            return mu1_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }
	
	//_AO_09/02/2014_
        const volScalarField& upsilon() const
        {
            return upsilon_;
        }
	
	//_AO_09/02/2014_
        const volSymmTensorField& tau() const
        //const volTensorField& tau() const
	{
            return tau_; // Convert symmTensor to Tensor, used in UEqns.H
	}

	//_AO_09/02/2014_
        const volScalarField& ppMagf() const
	{
            return ppMagf_; // Convert symmTensor to Tensor, used in UEqns.H
	}		
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
