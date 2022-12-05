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
    // Static data
    static const label pi_ = Foam::constant::mathematical::pi;

    static const Foam::dimensionedScalar Av_;
    static const Foam::dimensionedScalar kB_;
    static const Foam::dimensionedScalar Ru_;
    static const Foam::dimensionedScalar rho_soot_;
    static const Foam::dimensionedScalar W_carbon_;
    static const Foam::dimensionedScalar W_hydrogen_;
    static const Foam::dimensionedScalar C_min_;
    static const Foam::dimensionedScalar H_min_;
    static const Foam::dimensionedScalar PAH_rho_const_;

    // // Carbon and hydrogen atoms in the smallest soot particle
    // static const dimensioned<scalar> C_min_
    // (
    //     "C_min",
    //     dimensionSet(0,0,0,0,0,0,0),
    //     scalar(378)
    // );

    // static const dimensioned<scalar> H_min_
    // (
    //     "H_min_",
    //     dimensionSet(0,0,0,0,0,0,0),
    //     scalar(378)
    // );

    // // Constant for converting molecular weight to density
    // static const dimensioned<scalar> PAH_rho_const_
    // (
    //     "PAH_rho_const_",
    //     dimensionSet(0,-3,0,0,1,0,0),
    //     scalar(171943.5197)
    // );

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

    // Soot tracked fields
    volScalarField N_agg_;
    volScalarField N_pri_;
    volScalarField C_tot_;
    volScalarField H_tot_;

    // Morphology
    volScalarField n_p_;
    volScalarField d_p_;
    volScalarField d_m_;
    volScalarField d_g_;
    volScalarField A_tot_;

    // Species List
    wordList speciesList_;
    // Species Indicies
    HashTable<label> speciesIndicies_;

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

        virtual void findIndicies();

        virtual void updateMorphology();

        // Return Molecular Weight in kg/mol
        virtual dimensioned<scalar> W(label index)
        {
            basicSpecieMixture& composition = this->thermo().composition();
            return composition.W(index) / 1000 * dimensionedScalar("onekgPerMol", dimensionSet(1,0,0,0,-1,0,0), scalar(1));
        };
        
        // Return concentration in mol/m3
        virtual tmp<volScalarField> C(label index)
        {
            basicSpecieMixture& composition = this->thermo().composition();
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "C",
                    (composition.Y()[index] * this->thermo().rho() / W(index))
                )
            );
        };

        //- Return the number concentration of agglomerates N_agg
        // Units: mol/kg
        virtual tmp<volScalarField> N_agg() const
        {
            return N_agg_;
        };

        //- Return the number concentration of primary particles N_pri
        // Units: mol/kg
        virtual tmp<volScalarField> N_pri() const
        {
            return N_pri_;
        };

        //- Return the carbon contant of soot C_tot
        // Units: mol/kg
        virtual tmp<volScalarField> C_tot() const
        {
            return C_tot_;
        };

        //- Return the carbon contant of soot H_tot
        // Units: mol/kg
        virtual tmp<volScalarField> H_tot() const
        {
            return H_tot_;
        };

       
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
