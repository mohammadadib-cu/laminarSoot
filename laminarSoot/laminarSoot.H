/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::combustionModels::laminarSoot

Group
    grpCombustionModels

Description
    Laminar combustion model.

SourceFiles
    laminarSoot.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModels_laminarSoot_H
#define combustionModels_laminarSoot_H

#include "ChemistryCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class laminarSoot Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class laminarSoot
:
    public ChemistryCombustion<ReactionThermo>
{
    // Private data

        //- Integrate reaction rate over the time-step
        //  using the selected ODE solver
        bool integrateReactionRate_;

protected:

    // Protected Member Functions

        //- Return the chemical time scale
        tmp<volScalarField> tc() const;

    // Soot properties dictionary
    const IOdictionary sootProps_;

    // Inception Flag
    bool inception_enabled_;
    // Surface growth Flags
    bool HACA_growth_enabled_;
    bool PAH_growth_enabled_ ;
    bool use_alpha_emprical_;
    // Surface oxidation Flag
    bool oxidation_enabled_;
    // Coagulation flag
    bool coagulation_enabled_;

    // Precursor PAH Info
    wordList PAH_names_;
    labelList PAH_n_C_;
    labelList PAH_n_H_;
    labelList PAH_indicies_;

    // Dimer Info
    wordList dimer_names_;
    labelList dimer_n_C_;
    labelList dimer_n_H_;
    labelList dimer_PAH_1_index_;
    labelList dimer_PAH_2_index_;
    labelList dimer_PAH_1_id_;
    labelList dimer_PAH_2_id_;


private:

    // Private Member Functions

        //- No copy construct
        laminarSoot(const laminarSoot&) = delete;

        //- No copy assignment
        void operator=(const laminarSoot&) = delete;


public:

    //- Runtime type information
    TypeName("laminarSoot");


    // Constructors

        //- Construct from components
        laminarSoot
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~laminarSoot();


    // Member Data
    


    // Member Functions


        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();

        virtual bool readPAHs();

        virtual void buildDimers();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laminarSoot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
