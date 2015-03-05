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

Class
    Foam::radialModel

SourceFiles
    radialModel.C

\*---------------------------------------------------------------------------*/

#ifndef radialModel_H
#define radialModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class radialModel Declaration
\*---------------------------------------------------------------------------*/

class radialModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        radialModel(const radialModel&);

        //- Disallow default bitwise assignment
        void operator=(const radialModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("radialModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        radialModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        radialModel(const dictionary& dict);


    // Selectors

        static autoPtr<radialModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~radialModel();


    // Member Functions

        //- Radial distribution function
        virtual tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const = 0;

        //- Derivative of the radial distribution function
        virtual tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const = 0;

        //_AO_09/02/2014 - Radial distribution function with jamming effect (critical solid volume fraction)
        virtual tmp<volScalarField> g0jamming
        (
            const fvMesh& mesh,
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const scalar& alpha_f,
            const scalar& alpha_c
        ) const = 0;

        //_AO_09/08/2014 - Derivative of Radial distribution function with jamming effect 
        virtual tmp<volScalarField> g0jammingPrime
        (
            const fvMesh& mesh,
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const scalar& alpha_f,
            const scalar& alpha_c
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
