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
    Foam::kineticTheoryModels::radialModels::Torquato

Description

SourceFiles
    TorquatoRadial.C

\*---------------------------------------------------------------------------*/

#ifndef Torquato_H
#define Torquato_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{

/*---------------------------------------------------------------------------*\
                           Class Torquato Declaration
\*---------------------------------------------------------------------------*/

class Torquato
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("Torquato");


    // Constructors

        //- Construct from components
        Torquato(const dictionary& dict);


    //- Destructor
    virtual ~Torquato();


    // Member Functions


        tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;

        tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;

	//_AO_09/01/2014_Eq.12, p.4	
        tmp<volScalarField> g0jamming
        (
    	    const fvMesh& mesh,
	    const volScalarField& alpha,
    	    const dimensionedScalar& alphaMax,
            const scalar& alpha_f,
            const scalar& alpha_c
        ) const;

        tmp<volScalarField> g0jammingPrime
        (
    	    const fvMesh& mesh,
	    const volScalarField& alpha,
    	    const dimensionedScalar& alphaMax,
            const scalar& alpha_f,
            const scalar& alpha_c
        ) const;
			
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radialModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
