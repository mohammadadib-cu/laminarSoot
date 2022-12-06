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
    // HACA Surface growth Flags
    bool HACA_growth_enabled_;
    // Surface oxidation Flag
    bool HACA_oxidation_enabled_;
    // PAH Surface growth Flags
    bool PAH_growth_enabled_ ;
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

    // Source Terms
    // Inception
    volScalarField S_inc_N_;
    volScalarField S_inc_C_tot_;
    volScalarField S_inc_H_tot_;
    // Surface growth
    volScalarField S_grow_C_tot_;
    volScalarField S_grow_H_tot_;
    // Oxidation
    volScalarField S_ox_C_tot_;
    // Coagulation
    volScalarField S_coag_N_agg_;

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

        virtual void updateInception();

        virtual void updateGrowth();

        virtual void updateOxidation();
        
        virtual void updateCoagulation();

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
        //---------------------------------------------------------------------------------------------
        // Gas properties
        //- Returns mean free path of gas
        virtual tmp<volScalarField> lambda_gas() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    this->thermo().mu() / this->thermo().rho() * pow(pi_ * this->thermo().W() / (2.0* kB_ * Av_ * this->thermo().T()), 0.5)
                )
            );
        };

        //---------------------------------------------------------------------------------------------
        // Morphology
        //- Returns volume of each agglomerate 
        virtual tmp<volScalarField> V_agg()
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    C_tot() * W_carbon_ / (rho_soot_ * N_agg() * Av_)
                )
            );
        }; 
        //- Returns mass of each agglomerate 
        virtual tmp<volScalarField> m_agg()
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    C_tot() * W_carbon_ / (N_agg() * Av_)
                )
            );
        };   

        //---------------------------------------------------------------------------------------------
        // PAH Characteristics
        // Mass of PAH
        virtual dimensioned<scalar> m_PAH(label id)
        {
            label index = PAH_indicies_[id];
            return W(index) / Av_;            
        };

        // Diameter of PAH
        virtual dimensioned<scalar> d_PAH(label id)
        {
            const dimensionedScalar rho_PAH = PAH_rho_const_ * (W_carbon_ * PAH_n_C_[id] + W_hydrogen_ * PAH_n_H_[id]) / (PAH_n_C_[id] +  PAH_n_H_[id]);
            const dimensionedScalar V_PAH = m_PAH(id) / rho_PAH;
            return pow(6.0 * V_PAH / pi_, 1.0/3.0);
        };
        //---------------------------------------------------------------------------------------------
        // Dimerization
        // Return forward rate of physical dimerization in m3/mol-s
        virtual tmp<volScalarField> k_FWD(label id1, label id2)
        {
            const volScalarField& T = this->thermo().T();
            // reduced dimer diameter
            const dimensionedScalar d_reduced = 2 * d_PAH(id1) * d_PAH(id2) / (d_PAH(id1) + d_PAH(id2));
            // reduced dimer mass
            const dimensionedScalar m_reduced = m_PAH(id1) * m_PAH(id2) / (m_PAH(id1) + m_PAH(id2));

            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_FWD",
                    (1.0 * 2.2 * 0.1 * Av_ * d_reduced * d_reduced * pow(8.0 * pi_ * kB_ / m_reduced, 0.5) * pow(T, 0.5))
                )
            );
        };
        // Return reverse rate of physical dimerization in m3/mol-s
        virtual tmp<volScalarField> k_REV(label id1, label id2)
        {
            // Constants
            dimensionedScalar a_k_REV("a_k_REV", dimensionSet(0,0,0,0,0,0,0),  scalar(0.115));
            dimensionedScalar b_k_REV("b_k_REV", dimensionSet(0,0,0,0,0,0,0),  scalar(0.18));
            dimensionedScalar c_k_REV("c_k_REV", dimensionSet(0,2,-2,0,0,0,0), scalar(970688.0));
            dimensionedScalar d_k_REV("d_k_REV", dimensionSet(1,2,-2,0,-1,0,0),scalar(34069.0));

            const volScalarField& T = this->thermo().T();
            label index1 = PAH_indicies_[id1];
            label index2 = PAH_indicies_[id2];
            const dimensionedScalar W_reduced = W(index1) * W(index2) / (W(index1) + W(index2));
            const dimensionedScalar epsilon = c_k_REV * W_reduced - d_k_REV;

            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_REV",
                    (k_FWD(id1, id2) * pow(10.0, -b_k_REV) * exp(-a_k_REV * epsilon * 2.3025851/(Ru_ * T)))
                )
            );
        };

        // Return rate of chemical dimerization in m3/mol-s
        virtual tmp<volScalarField> k_REAC()
        {
            // Chemical dimerization rate constants
            dimensionedScalar Ea_k_REAC("Ea_k_REAC", dimensionSet(1,2,-2,0,-1,0,0), scalar(96232.0));
            dimensionedScalar A_k_REAC("A_k_REAC", dimensionSet(0,3,-1,0,-1,0,0), scalar(5.0e6));
            const volScalarField& T = this->thermo().T();
            
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_REAC",
                    (A_k_REAC * exp(-Ea_k_REAC/(Ru_ * T)))
                )
            );
        }

        // Returns dimer formation mol/m3-s
        virtual tmp<volScalarField> dimerROP(label id1, label id2)
        {
            label index1 = PAH_indicies_[id1];
            label index2 = PAH_indicies_[id2];
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "dimerROP",
                    (k_REAC() * k_FWD(id1, id2) * C(index1) * C(index2) / (k_REAC() + k_REV(id1, id2)))
                )
            );
        }

        //---------------------------------------------------------------------------------------------
        // HACA Rate
        virtual tmp<volScalarField> k_4_HACA()
        {
            const volScalarField& T = this->thermo().T();
            dimensionedScalar oneKelvin(dimTemperature, scalar(1.0));
            dimensionedScalar k_unit(dimensionSet(0,3,-1,0,-1,0,0), scalar(1.0));
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_4_HACA",
                    (k_unit * 8.00e1 * pow(T/oneKelvin, 1.56) * exp(oneKelvin * (-1912.43) / T))
                )
            );
        }
        virtual tmp<volScalarField> k_5_HACA()
        {
            const volScalarField& T = this->thermo().T();
            dimensionedScalar oneKelvin(dimTemperature, scalar(1.0));
            dimensionedScalar k_unit(dimensionSet(0,3,-1,0,-1,0,0), scalar(1.0));
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_5_HACA",
                    (k_unit * 2.20e6 * exp(oneKelvin * (-3774.53) / T))
                )
            );
        } 
        virtual tmp<volScalarField> k_6_HACA()
        {
            const volScalarField& T = this->thermo().T();
            dimensionedScalar oneKelvin(dimTemperature, scalar(1.0));
            dimensionedScalar k_unit(dimensionSet(0,3,-1,0,-1,0,0), scalar(1.0));
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "k_6_HACA",
                    (k_unit * 0.13e0 * pow(T/oneKelvin, 0.0))
                )
            );
        }       
        virtual tmp<volScalarField> C_soot_0()
        {
            // Species indicies
            label H_i = speciesIndicies_["H"];
            label OH_i = speciesIndicies_["OH"];
            label H2_i = speciesIndicies_["H2"];
            label H2O_i = speciesIndicies_["H2O"];
            label C2H2_i = speciesIndicies_["C2H2"];
            label O2_i = speciesIndicies_["O2"];
            
            const volScalarField& T = this->thermo().T();

            dimensionedScalar oneKelvin(dimTemperature, scalar(1.0));
            dimensionedScalar k_unit(dimensionSet(0,3,-1,0,-1,0,0), scalar(1.0));

            // Forward Reaction Rates
            // Units : m3/mol-s
            volScalarField k_1(k_unit * 4.17e7 * exp(oneKelvin * (-6542.52) / T));
            volScalarField k_2(k_unit * 1.00e4 * pow(T/oneKelvin, 0.734) * exp(oneKelvin * (-719.68) / T));
            volScalarField k_3(k_unit * 2.00e7 * pow(T/oneKelvin, 0.0));
            // volScalarField k_4(k_unit * 8.00e1 * pow(T/oneKelvin, 1.56) * exp(oneKelvin * (-1912.43) / T));
            // volScalarField k_5(k_unit * 2.20e6 * exp(oneKelvin * (-3774.53) / T));
            // volScalarField k_6(k_unit * 0.13e0 * pow(T/oneKelvin, 0.0));

            // Reverse Reaction Rates
            // Units : m3/mol-s
            volScalarField k_r_1(k_unit * 3.9e6 * exp(oneKelvin * (-5535.98) / T));
            volScalarField k_r_2(k_unit * 3.68e2 * pow(T/oneKelvin, 1.139) * exp(oneKelvin * (-8605.94) / T));

            // chi_soot_0
            dimensionedScalar chi_soot_0_denum_limit("chi_soot_0_denum_limit",dimless/dimTemperature, scalar(1e-20));
            volScalarField chi_soot_0_denum
            (
                max 
                (
                    k_r_1 * C(H2_i) + k_r_2 * C(H2O_i) + k_3 * C(H_i) + k_4_HACA() * C(C2H2_i) + k_5_HACA() * C(O2_i) + k_1 * C(H_i) + k_2 * C(OH_i),
                    chi_soot_0_denum_limit
                )
            );

            dimensionedScalar chi_soot_CH(dimensionSet(0,-2,0,0,0,0,0), scalar(2.3e19));
            volScalarField chi_soot_0
            (
                (k_1 * C(H_i) + k_2 * C(OH_i)) / chi_soot_0_denum * chi_soot_CH   
            );

            // C_soot_0
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "C_soot_0",
                    (A_tot_ / Av_ * chi_soot_0)
                )
            );
        }

        virtual tmp<volScalarField> alpha()
        {
            dimensionedScalar oneKelvin(dimTemperature, scalar(1.0));
            dimensionedScalar alpha_min(dimless, scalar(0.0));
            dimensionedScalar alpha_max(dimless, scalar(1.0));
            const volScalarField& T = this->thermo().T();
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    min
                    (
                        max
                        (
                            tanh
                            (
                                (12.56 - 0.00563 * T / oneKelvin) / log10 ( rho_soot_ * pi_ / 6.0 * pow(d_p_, 3.0) * Av_ / W_carbon_ ) - 
                                1.38 + 0.00068 * T / oneKelvin
                            ),
                            alpha_min
                        ),
                        alpha_max
                    )
                )
            );
        }

        virtual tmp<volScalarField> HACAGrowthRate()
        {
            label C2H2_i = speciesIndicies_["C2H2"];
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    2 * alpha() * k_4_HACA() * C(C2H2_i) * C_soot_0()
                )
            );
        }

        virtual tmp<volScalarField> HACAOxidationRate()
        {
            label OH_i = speciesIndicies_["OH"];
            label O2_i = speciesIndicies_["O2"];
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    -2 * alpha() * k_5_HACA() * C(O2_i) * C_soot_0() - k_6_HACA() * C(OH_i) * N_agg()
                )
            );
        }


        //---------------------------------------------------------------------------------------------
        // PAH Adsorption Rate

        virtual tmp<volScalarField> PAHAdsorptionRate(label id)
        {

            // PAH Index and Id
            label index = PAH_indicies_[id];
            // Temperature
            const volScalarField& T = this->thermo().T();
            // Viscosity of gas
            volScalarField mu = this->thermo().mu();
            // beta fm coag PAH-soot
            volScalarField beta_fm_soot_PAH
            (
                2.2 * pow(pi_ * kB_ * T / 2.0 * (1.0/m_agg() + 1.0/m_PAH(id)), 0.5) * pow(d_g_ + d_PAH(id), 2.0)
            );
            // beta cont coag PAH-soot
            volScalarField C_s_soot_d_m
            (
                1.0 + (2.0 * lambda_gas() / d_m_ ) * (1.21 + 0.4*exp(-0.78 * d_m_ / lambda_gas()))
            );
            volScalarField C_s_soot_d_PAH
            (
                1.0 + (2.0 * lambda_gas() / d_PAH(id)) * (1.21 + 0.4*exp(-0.78 * d_PAH(id) / lambda_gas()))
            );
            volScalarField beta_cont_soot_PAH
            (
                2 * kB_ * T / (3 * mu) * (C_s_soot_d_m / d_g_ + C_s_soot_d_PAH / d_PAH(id)) * (d_g_ + d_PAH(id))
            );
            // Forward reaction rate
            volScalarField k_FWD_soot_PAH
            (
                beta_fm_soot_PAH * beta_cont_soot_PAH / (beta_fm_soot_PAH + beta_cont_soot_PAH) * Av_
            );
            // W reduced
            // Constants
            dimensionedScalar a_k_REV_soot_PAH(dimensionSet(0,0,0,0,0,0,0),  scalar(0.115));
            dimensionedScalar b_k_REV_soot_PAH(dimensionSet(0,0,0,0,0,0,0),  scalar(0.18));
            dimensionedScalar c_k_REV_soot_PAH(dimensionSet(0,2,-2,0,0,0,0), scalar(970688.0));
            dimensionedScalar d_k_REV_soot_PAH(dimensionSet(1,2,-2,0,-1,0,0),scalar(34069.0));
            // W soot
            volScalarField W_soot(C_tot() * W_carbon_ / N_agg());
            volScalarField epsilon_soot_PAH
            (
                c_k_REV_soot_PAH * (W(index) * W_soot) / (W(index) + W_soot) - d_k_REV_soot_PAH
            );
            // The reverse rate of physical dimerization
            volScalarField k_REV_soot_PAH
            (
                k_FWD_soot_PAH * pow(10.0, -b_k_REV_soot_PAH) * exp(-a_k_REV_soot_PAH * epsilon_soot_PAH * 2.3025851/(Ru_ * T))
            );

            // Chemical adsorption rate
            dimensionedScalar Ea_k_REAC_soot_PAH(dimensionSet(1,2,-2,0,-1,0,0), scalar(209631.0));
            dimensionedScalar A_k_REAC_soot_PAH(dimensionSet(0,3,-1,0,-1,0,0), scalar(1.0e10));
            volScalarField k_REAC_soot_PAH
            (
                A_k_REAC_soot_PAH * exp(-Ea_k_REAC_soot_PAH / (Ru_ * T))
            );
            // Soot concentration
            // volScalarField C_soot(N_agg * rho);
            // k_REAC_soot_PAH * k_FWD_soot_PAH * C(index) * (N_agg * rho) / (k_REAC_soot_PAH + k_REV_soot_PAH)
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    k_REAC_soot_PAH * k_FWD_soot_PAH * C(index) * (N_agg() * this->thermo().rho()) / (k_REAC_soot_PAH + k_REV_soot_PAH)
                )
            );
        }

        //---------------------------------------------------------------------------------------------
        // Coagulation

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
