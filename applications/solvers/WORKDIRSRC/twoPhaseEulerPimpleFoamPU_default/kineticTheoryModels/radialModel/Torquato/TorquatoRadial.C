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

#include "TorquatoRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(Torquato, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        Torquato,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Torquato::Torquato
(
    const dictionary& dict
)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::Torquato::~Torquato()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Torquato::g0
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
Foam::kineticTheoryModels::radialModels::Torquato::g0prime
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
Foam::kineticTheoryModels::radialModels::Torquato::g0jamming
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const scalar& alpha_f,
    const scalar& alpha_c
) const
{
    
    const scalar constSMALL = 1.e-06;
    volScalarField valg0(alpha);
    
    forAll(mesh.cells(),ii)
    {
	if( alpha[ii] <= alpha_f )
	{
	   valg0[ii] =   (1.0 - alpha[ii]/2.) / (pow(1.0 - alpha[ii], 3));
	}else
	{
	   valg0[ii] =   (1.0 - alpha[ii]/2.) / (pow(1.0 - alpha[ii], 3))
	               * (alpha_c - alpha_f)
		       / (alpha_c - alpha[ii] + constSMALL ) ;    
	}
    }

    return  valg0;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::Torquato::g0jammingPrime
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const scalar& alpha_f,
    const scalar& alpha_c
) const
{
    return
        2.5/sqr(1.0 - alpha)
      + 4.0*alpha/pow(1.0 - alpha, 3.0)
      + 1.5*sqr(alpha)/pow(1.0 - alpha, 4.0);
}

// ************************************************************************* //
