/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    defineTypeNameAndDebug(ChialvoSundaresanRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        ChialvoSundaresanRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanRadial::ChialvoSundaresanRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanRadial::~ChialvoSundaresanRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0
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


Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return 
        - alpha/sqr(1.0 - alpha) 
        + (3.0*(1.0 - alpha) + 6.0*sqr(alpha))/(2.0*(1.0 - alpha))
        + (2.0*alpha*(1.0 - alpha) + 3.0*pow(alpha, 3))
         /(2.0*pow(1.0 - alpha, 4));
}

//_AO_09/08/2014
Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0jamming
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
    
    volScalarField valg0CS = (1.0 - alpha/2.) / (pow(1.0 - alpha, 3));
    
    return  
    	valg0CS +  alpha2 * alpha * alpha 
	         / pow( alpha_c - alpha , 1.5);
}

Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0jammingPrime
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
    
    volScalarField valg0CSprime = 
    - 1.0 / ( 2.0 * pow(1.0 - alpha, 3))
    - 3.0 * ( 1.0 - alpha/2.) / (pow(1.0 - alpha, 4));
    
    return  
    	valg0CSprime  
    + 2.0 * alpha2 * alpha / pow( alpha_c - alpha ,1.5)
    - 3.0 / 2.0 * alpha2 * alpha * alpha / pow( alpha_c - alpha, 2.5);

}


// ************************************************************************* //
