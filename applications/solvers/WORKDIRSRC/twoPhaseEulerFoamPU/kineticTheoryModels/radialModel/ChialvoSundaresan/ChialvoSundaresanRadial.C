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

\*---------------------------------------------------------------------------*/

#include "ChialvoSundaresanRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(ChialvoSundaresan, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        ChialvoSundaresan,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::ChialvoSundaresan
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::~ChialvoSundaresan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    return
        1.0/(1.0 - alpha)
      + 3.0*alpha/(2.0*sqr(1.0 - alpha))
      + sqr(alpha)/(2.0*pow(1.0 - alpha, 3));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return
        2.5/sqr(1.0 - alpha)
      + 4.0*alpha/pow(1.0 - alpha, 3.0)
      + 1.5*sqr(alpha)/pow(1.0 - alpha, 4.0);
}

//_AO_09/01/2014_Eq.12, p.4
Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0jamming
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const scalar& alpha_f,
    const scalar& alpha_c
) const
{
    //AO_09/02/2014 Eq. 31, p. 10
    const scalar alpha2 = 0.58;
    const scalar scalarSMALL = 1.e-06;
    
    volScalarField valg0CS = (1.0 - alpha/2.) / (pow(1.0 - alpha, 3));
    
    return  
    	valg0CS +  alpha2 * alpha * alpha 
	         / pow((alpha_c - alpha + scalarSMALL),1.5);
	
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::ChialvoSundaresan::g0jammingPrime
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const scalar& alpha_f,
    const scalar& alpha_c
) const
{
    //AO_09/02/2014 Eq. 31, p. 10
    const scalar alpha2 = 0.58;
    const scalar scalarSMALL = 1.e-06;
    
    volScalarField valg0CSprime = 
    - 1.0 / ( 2.0 * pow(1.0 - alpha, 3))
    - 3.0 * ( 1.0 - alpha/2.) / (pow(1.0 - alpha, 4));
    
    return  
    	valg0CSprime  
    + 2.0 * alpha2 * alpha / pow((alpha_c - alpha + scalarSMALL),1.5)
    - 3.0 / 2.0 * alpha2 * alpha * alpha / pow((alpha_c - alpha + scalarSMALL), 2.5);
	
}

// ************************************************************************* //
