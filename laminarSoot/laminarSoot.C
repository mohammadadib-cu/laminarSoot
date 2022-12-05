/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "laminarSoot.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"
#include "IFstream.H"


// Static data
template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::Av_ =     
    Foam::dimensionedScalar(
        "Av",
        Foam::dimensionSet(0,0,0,0,-1,0,0),
        scalar(6.0221409e+23)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::kB_ =     
    Foam::dimensionedScalar(
        "kB",
        Foam::dimensionSet(1,2,-2,-1,0,0,0),
        scalar(1.38064852e-23)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::Ru_ =     
    Foam::dimensionedScalar(
        "Ru_",
        Foam::dimensionSet(1,2,-2,-1,-1,0,0),
        scalar(8.314462618)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::rho_soot_ =     
    Foam::dimensionedScalar(
        "rho_soot_",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        scalar(1800.0)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::W_carbon_ =     
    Foam::dimensionedScalar(
        "W_carbon",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        scalar(0.012)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::W_hydrogen_ =     
    Foam::dimensionedScalar(
        "W_hydrogen_",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        scalar(0.001)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::C_min_ =     
    Foam::dimensionedScalar(
        "C_min",
        Foam::dimensionSet(0,0,0,0,0,0,0),
        scalar(378)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::H_min_ =     
    Foam::dimensionedScalar(
        "H_min_",
        Foam::dimensionSet(0,0,0,0,0,0,0),
        scalar(378)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::PAH_rho_const_ =     
    Foam::dimensionedScalar(
        "PAH_rho_const_",
        Foam::dimensionSet(0,-3,0,0,1,0,0),
        scalar(171943.5197)
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



template<class ReactionThermo>
Foam::combustionModels::laminarSoot<ReactionThermo>::laminarSoot
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ChemistryCombustion<ReactionThermo>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    integrateReactionRate_
    (
        this->coeffs().getOrDefault("integrateReactionRate", true)
    ),
    sootProps_
    (
        IOobject
        (
            "sootProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    inception_enabled_(true),
    HACA_growth_enabled_(true),
    PAH_growth_enabled_(true),
    use_alpha_emprical_(true),
    oxidation_enabled_(true),
    coagulation_enabled_(true),
    N_agg_
    (
        IOobject
        (
            "N_agg",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    N_pri_
    (
        IOobject
        (
            "N_pri",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    C_tot_
    (
        IOobject
        (
            "C_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    H_tot_
    (
        IOobject
        (
            "H_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh()
    ),
    n_p_
    (
        IOobject
        (
            "n_p",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),dimensionedScalar("n_p", dimensionSet(0,0,0,0,0,0,0),1.0)
    ),
    d_p_
    (
        IOobject
        (
            "d_p",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),dimensionedScalar("d_p", dimensionSet(0,1,0,0,0,0,0), 2.0e-9)
    ),
    d_m_
    (
        IOobject
        (
            "d_m",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),dimensionedScalar("d_m", dimensionSet(0,1,0,0,0,0,0), 2.0e-9)
    ),
    d_g_
    (
        IOobject
        (
            "d_g",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),dimensionedScalar("d_g", dimensionSet(0,1,0,0,0,0,0), 2.0e-9)
    ),
    A_tot_
    (
        IOobject
        (
            "A_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),dimensionedScalar("A_tot", dimensionSet(-1,2,0,0,0,0,0), 2.0e-9)
    )

{

    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }

    speciesList_ = {"H", "H2", "OH", "H2O", "C2H2", "O2", "CO"};
    findIndicies();
    readPAHs();
    buildDimers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminarSoot<ReactionThermo>::~laminarSoot()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminarSoot<ReactionThermo>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::correct()
{
    if (this->active())
    {
        if (integrateReactionRate_)
        {
            if (fv::localEulerDdt::enabled(this->mesh()))
            {
                const scalarField& rDeltaT =
                    fv::localEulerDdt::localRDeltaT(this->mesh());

                scalar maxTime;
                if (this->coeffs().readIfPresent("maxIntegrationTime", maxTime))
                {
                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
        updateMorphology();
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminarSoot<ReactionThermo>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    fvScalarMatrix& Su = tSu.ref();

    if (this->active())
    {
        const label specieI =
            this->thermo().composition().species()[Y.member()];

        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminarSoot<ReactionThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                this->thermo().phasePropertyName(typeName + ":Qdot"),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    if (this->active())
    {
        tQdot.ref() = this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class ReactionThermo>
bool Foam::combustionModels::laminarSoot<ReactionThermo>::read()
{
    if (ChemistryCombustion<ReactionThermo>::read())
    {
        integrateReactionRate_ =
            this->coeffs().getOrDefault("integrateReactionRate", true);
        return true;
    }

    

    return false;
}

template<class ReactionThermo>
bool Foam::combustionModels::laminarSoot<ReactionThermo>::readPAHs()
{
    sootProps_.readEntry("PAHs", PAH_names_);
    PAH_n_C_.resize(PAH_names_.size());
    PAH_n_H_.resize(PAH_names_.size());
    PAH_indicies_.resize(PAH_names_.size());

    const ReactionThermo& thermo = this->thermo();
    const dictionary thermoDict = IFstream(fileName(thermo.lookup("foamChemistryThermoFile")).expand())();

    forAll(PAH_names_, i)
    {
        PAH_indicies_[i] = this->thermo().composition().species()[PAH_names_[i]];
        const dictionary* elemsDict = thermoDict.subDict(PAH_names_[i]).findDict("elements");
        wordList elemNames(elemsDict->toc());
        
        forAll(elemNames, eni)
        {
            if (elemNames[eni] == "C")
            {
                PAH_n_C_[i] = elemsDict->getOrDefault<label>
                (
                    elemNames[eni],
                    0
                );        	
            }else if (elemNames[eni] == "H"){
                PAH_n_H_[i] = elemsDict->getOrDefault<label>
                (
                    elemNames[eni],
                    0
                );         	
            
            }

        }

    };

    // Outputing the list of PAHs
    Info << "Incepient species\n" << endl;
    forAll(PAH_names_, i)
    {
        const label index = PAH_indicies_[i];
        Info << PAH_names_[i] << ": C" << PAH_n_C_[i] << "H" << PAH_n_H_[i] << " index in species: " << PAH_indicies_[i] << " MW:" << this->thermo().composition().W(index) << endl;
    }

    return true;
}

// Building Dimers
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::buildDimers()
{
    basicSpecieMixture& composition = this->thermo().composition();

    label dimer_size = PAH_names_.size() * (PAH_names_.size() +1 ) / 2;
    Info << dimer_size << " possible dimers are created! \n" << endl;
    dimer_names_.resize(dimer_size);
    dimer_n_C_.resize(dimer_size);
    dimer_n_H_.resize(dimer_size);
    dimer_PAH_1_index_.resize(dimer_size);
    dimer_PAH_2_index_.resize(dimer_size);
    dimer_PAH_1_id_.resize(dimer_size);
    dimer_PAH_2_id_.resize(dimer_size);


    label dimer_id = -1;

    for (int i = 0; i < PAH_names_.size(); i++) {
        for (int j = i; j < PAH_names_.size(); j++) {
            dimer_id += 1;
            dimer_names_[dimer_id] = PAH_names_[i] + PAH_names_[j];
            dimer_n_C_[dimer_id] = PAH_n_C_[i] + PAH_n_C_[j];
            dimer_n_H_[dimer_id] = PAH_n_H_[i] + PAH_n_H_[j];
            dimer_PAH_1_index_[dimer_id] = composition.species()[PAH_names_[i]];
            dimer_PAH_2_index_[dimer_id] = composition.species()[PAH_names_[j]];
            dimer_PAH_1_id_[dimer_id] = i;
            dimer_PAH_2_id_[dimer_id] = j;
        }
    }



    Info << "Outputing Dimers \n" << endl;
    forAll(dimer_names_, i)
    {
        Info << dimer_names_[i] << ": C" <<  dimer_n_C_[i] << "H" << dimer_n_H_[i] 
        << " PAH1 index:" << dimer_PAH_1_index_[i] << " PAH2 index:" << dimer_PAH_2_index_[i] << endl;
    }

}

// Building Dimers

template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::findIndicies()
{
    basicSpecieMixture& composition = this->thermo().composition();
    forAll(speciesList_, i)
    {
        const label specieIndex = composition.species()[speciesList_[i]];
        if (!composition.species().found(speciesList_[i]))
        {
            FatalIOErrorIn("laminarSoot::findIndicies()", this->thermo())
                << speciesList_[i] <<" is not found in available species "
                << composition.species() << exit(FatalIOError);
        }
        speciesIndicies_.insert
        (
            speciesList_[i],
            specieIndex
        );
        Info << speciesList_[i] << " is found!" << endl;
    }
}


// Updating morphology
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateMorphology()
{
    // number of primary particles
    n_p_ = N_pri() / N_agg();
    n_p_.max(1.0);
    if (!coagulation_enabled_)
    {
        Info<< "enforcing n_p = 1\n" << endl;
        n_p_.min(1.00000001);
    }

    // Total volume of soot per kg of gas
    volScalarField V_tot(C_tot() * W_carbon_ / rho_soot_);
    // Volume of each primary particles
    volScalarField V_p(V_tot / (N_pri() * Av_));
    // Volume of each agglomerate
    volScalarField V_agg(V_tot / (N_agg() * Av_));
    // Mass of agglomerate
    volScalarField m_agg(V_agg * rho_soot_); 

    // Primary particle diameter
    d_p_ = pow(6.0 * V_p / (pi_), 1.0/3.0);

    // Surface area of each primary particle
    volScalarField A_p(pi_ * d_p_ * d_p_);

    // Total surface area
    A_tot_ = N_pri() * Av_ * A_p;

    //  Mobility diameter
    d_m_ = d_p_ * pow(n_p_, 0.45);

    // Gyration Diameter
    volScalarField n_p_lowerlimit (n_p_*0.0+1.5);
    d_g_ = (n_p_ <= n_p_lowerlimit) * (d_m_ / 1.29) + (n_p_ > n_p_lowerlimit) * (d_m_ / (pow(n_p_, -0.2)+0.4));

}


// ************************************************************************* //
