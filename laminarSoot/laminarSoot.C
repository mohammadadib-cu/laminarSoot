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
#include "fvOptions.H"


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
        scalar(12.011e-3)
    );

template<class ReactionThermo> const Foam::dimensionedScalar 
Foam::combustionModels::laminarSoot<ReactionThermo>::W_hydrogen_ =     
    Foam::dimensionedScalar(
        "W_hydrogen_",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        scalar(1.00784e-3)
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
        scalar(20)
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
    inception_enabled_(sootProps_.getOrDefault("inception_enabled", true)),
    HACA_growth_enabled_(sootProps_.getOrDefault("HACA_growth_enabled", true)),
    HACA_oxidation_enabled_(sootProps_.getOrDefault("HACA_oxidation_enabled", true)),
    PAH_growth_enabled_(sootProps_.getOrDefault("PAH_growth_enabled", true)),
    coagulation_enabled_(sootProps_.getOrDefault("coagulation_enabled", true)),
    scrubbing_enabled_(sootProps_.getOrDefault("scrubbing_enabled", true)),
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
    ),
    S_inc_N_
    (
        IOobject
        (
            "S_inc_N",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_inc_N", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_inc_C_tot_
    (
        IOobject
        (
            "S_inc_C_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_inc_C_tot", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_inc_H_tot_
    (
        IOobject
        (
            "S_inc_H_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_inc_H_tot", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_grow_C_tot_
    (
        IOobject
        (
            "S_grow_C_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_grow_C_tot", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_grow_H_tot_
    (
        IOobject
        (
            "S_grow_H_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_grow_H_tot", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_ox_C_tot_(
        IOobject
        (
            "S_ox_C_tot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_ox_C_tot", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    S_coag_N_agg_
    (
        IOobject
        (
            "S_coag_N_agg",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(), dimensionedScalar("S_coag_N_agg", dimensionSet(-1,0,-1,0,1,0,0),1.0e-30 )
    ),
    PAH_names_
    (
        sootProps_.lookup("PAHs")
    ),
    HACASpeciesList_(
        {"H", "H2", "OH", "H2O", "C2H2", "O2", "CO"}
    ),
    speciesList_
    (
        PAH_names_.size() + HACASpeciesList_.size()
    ),
    SR_(
        speciesList_.size()
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


    createPAHProps();
    createDimerProps();
    createSpeciesProps();
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
        resetSR();
        updateMorphology();
        updateInception();
        updateGrowth();
        updateOxidation();
        updateCoagulation();
        updateSoot();
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

        if (scrubbing_enabled_){
            if (speciesList_.found(Y.member()))
            {
                label spid = speciesIds_[Y.member()];
                Su += SR_[spid];
            }
        }
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
        if (scrubbing_enabled_){

            scalarField& Qdot = tQdot.ref();
            
            // basicSpecieMixture& composition = this->thermo().composition();
            // PtrList<volScalarField>& Y = this->thermo().composition().Y();
            forAll(this->thermo().composition().Y(), i)
            {
                if (speciesList_.found(this->thermo().composition().Y()[i].member()))
                {
                    label spid = speciesIds_[this->thermo().composition().Y()[i].member()];
                    forAll(Qdot, celli)
                    {
                        const scalar hi = this->thermo().composition().Hc(i);
                        Qdot[celli] -= hi*(this->chemistryPtr_->RR(i)[celli]+SR_[spid][celli]);
                    }
                }else{
                    forAll(Qdot, celli)
                    {
                        const scalar hi = this->thermo().composition().Hc(i);
                        Qdot[celli] -= hi*this->chemistryPtr_->RR(i)[celli];
                    }
                }
            }
        }else{
            tQdot.ref() = this->chemistryPtr_->Qdot();
        } 
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
bool Foam::combustionModels::laminarSoot<ReactionThermo>::createPAHProps()
{
    Info << "PAH names: " << PAH_names_ << endl;
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
void Foam::combustionModels::laminarSoot<ReactionThermo>::createDimerProps()
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
void Foam::combustionModels::laminarSoot<ReactionThermo>::createSpeciesProps()
{
    basicSpecieMixture& composition = this->thermo().composition();
    
    forAll(PAH_names_,i){
        speciesList_[i] = PAH_names_[i];
    }
    forAll(HACASpeciesList_,i){
        label spid = i + PAH_names_.size();
        speciesList_[spid]= HACASpeciesList_[i];
    }
    Info << "List of species: " << speciesList_ << endl;


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
        speciesIds_.insert
        (
            speciesList_[i],
            i
        );
        Info << speciesList_[i] << " is found! Index= " << speciesIndicies_[speciesList_[i]] << " Id: " << speciesIds_(speciesList_[i]) << endl;
    }
    
    // Create the fields for the chemistry sources
    forAll(SR_, fieldi)
    {
        SR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "SR." + speciesList_[fieldi],
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
            )
        );
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
    // Primary particle diameter
    d_p_ = pow (
        (6.0 / pi_) *
        (C_tot() * W_carbon_) / rho_soot_ *
        1.0 / (N_pri() * Av_)
        , 1.0/3.0
    );
    // Surface area of each primary particle
    A_tot_ = N_pri() * Av_ * pi_ * d_p_ * d_p_;
    //  Mobility diameter
    d_m_ = d_p_ * pow(n_p_, 0.45);
    // Gyration Diameter
    volScalarField n_p_lowerlimit (n_p_*0.0+1.5);
    d_g_ = (n_p_ <= n_p_lowerlimit) * (d_m_ / 1.29) + (n_p_ > n_p_lowerlimit) * (d_m_ / (pow(n_p_, -0.2)+0.4));

}

// Updating inception source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateInception()
{
    S_inc_N_ *= 0.0;
    S_inc_C_tot_ *= 0.0;
    S_inc_H_tot_ *= 0.0;
    if (inception_enabled_){
        volScalarField rho = this->thermo().rho();
        forAll(dimer_names_, i)
        {
            // PAH Index and Id
            label id1 = dimer_PAH_1_id_[i];
            label id2 = dimer_PAH_2_id_[i];
            volScalarField dimerROPField(dimerROP(id1, id2));
            // N_agg Source Term
            S_inc_N_ += (dimer_n_C_[i] / C_min_) * dimerROPField / rho;
            // C_tot Source Term
            S_inc_C_tot_ += dimer_n_C_[i] * dimerROPField / rho;
            // H_tot Source Term
            S_inc_H_tot_ += dimer_n_H_[i] * dimerROPField / rho;

            if (scrubbing_enabled_){
                // PAHs
                // species id
                label spid1 = speciesIds_[PAH_names_[id1]];
                label spid2 = speciesIds_[PAH_names_[id2]];
                // species index
                label spindex1 = speciesIndicies_[PAH_names_[id1]];
                label spindex2 = speciesIndicies_[PAH_names_[id2]];       
                SR_[spid1] -= dimerROPField * W(spindex1);
                SR_[spid2] -= dimerROPField * W(spindex2);

                // H2
                label H2_id = speciesIds_["H2"];
                label H2_index = speciesIndicies_["H2"];     
                SR_[H2_id] += dimerROPField * W(H2_index);
            }
        }
    }
}

// Updating growth source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateGrowth()
{
    S_grow_C_tot_ *= 0.0;
    S_grow_H_tot_ *= 0.0;

    if (HACA_growth_enabled_)
    {
        volScalarField rho = this->thermo().rho();
        volScalarField HACAGrowthRateField(HACAGrowthRate());
        S_grow_C_tot_ += 2 * HACAGrowthRateField / rho;
        S_grow_H_tot_ += 2 * HACAGrowthRateField * (0.25 / 2.00) / rho;

        if (scrubbing_enabled_){
            // C2H2
            label C2H2_id = speciesIds_["C2H2"];
            label C2H2_index = speciesIndicies_["C2H2"];             
            SR_[C2H2_id] -= HACAGrowthRateField * W(C2H2_index);

            // H
            label H_id = speciesIds_["H"];
            label H_index = speciesIndicies_["H"];             
            SR_[H_id] += HACAGrowthRateField * W(H_index) * (1.75 / 2.00);
        }
    }
    if (PAH_growth_enabled_)
    {
        volScalarField rho = this->thermo().rho();
        forAll(PAH_names_, id)
        {
            volScalarField PAHAdsorptionRateField(PAHAdsorptionRate(id));
            S_grow_C_tot_ += PAH_n_C_[id] * PAHAdsorptionRateField / rho;
            S_grow_H_tot_ += (PAH_n_H_[id] - 2) * PAHAdsorptionRateField / rho;

            if (scrubbing_enabled_){
                // PAH
                // species id
                label spid = speciesIds_[PAH_names_[id]];
                // species index
                label spindex = speciesIndicies_[PAH_names_[id]];   
                SR_[spid] -= PAHAdsorptionRateField * W(spindex);

                // H
                label H_id = speciesIds_["H"];
                label H_index = speciesIndicies_["H"];             
                SR_[H_id] += PAHAdsorptionRateField * W(H_index) * 2;
            }
        }
    }

}

// Updating oxidation source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateOxidation()
{
    S_ox_C_tot_ *= 0.0;
    if (HACA_oxidation_enabled_){
        volScalarField rho(this->thermo().rho());
        volScalarField HACAO2OxidationRateField(HACAO2OxidationRate());
        volScalarField HACAOHOxidationRateField(HACAOHOxidationRate());
        S_ox_C_tot_ += -1 * (HACAO2OxidationRateField + HACAOHOxidationRateField) / rho;

        if (scrubbing_enabled_){
            // O2
            label O2_id = speciesIds_["O2"];
            label O2_index = speciesIndicies_["O2"];     
            SR_[O2_id] -= 0.5 * HACAO2OxidationRateField * W(O2_index);

            // CO2
            label CO_id = speciesIds_["CO"];
            label CO_index = speciesIndicies_["CO"];     
            SR_[CO_id] += (HACAO2OxidationRateField + HACAOHOxidationRateField) * W(CO_index);

            // OH
            label OH_id = speciesIds_["OH"];
            label OH_index = speciesIndicies_["OH"];     
            SR_[OH_id] -= HACAOHOxidationRateField * W(OH_index);
        }
    }
}

// Updating coagulation source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateCoagulation()
{
    const volScalarField& T = this->thermo().T();
    volScalarField mu (this->thermo().mu());
    volScalarField rho (this->thermo().rho());
    // Free Molecule
    volScalarField beta_fm
    (
        4 * pow(pi_ * kB_ * T / m_agg(), 0.5) * d_g_ * d_g_
    );
    // Continuum
    volScalarField beta_cont(
        (8 * kB_ / (3 * mu)) * T * ( 1.0 + (2.0 * lambda_gas() / d_m_ )*(1.21 + 0.4*exp(-0.78*d_m_/lambda_gas())))
    );
    // Coagulation source term
    if (coagulation_enabled_){
        S_coag_N_agg_ = 0.5 * 1.82 * beta_fm * beta_cont / (beta_fm + beta_cont) * pow(N_agg(), 2.0) * Av_ * rho;   
    }
}

// Updating coagulation source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::resetSR()
{
    if (scrubbing_enabled_){
        forAll(SR_, fieldi){
            forAll(SR_[fieldi], celli){
                SR_[fieldi][celli] = 0.0;
            }
        }
    }
}

// Updating coagulation source terms
template<class ReactionThermo>
void Foam::combustionModels::laminarSoot<ReactionThermo>::updateSoot()
{
    const surfaceScalarField& phi = this->phi();
    const volScalarField D(diffusionCoeff());
    volScalarField rho (this->thermo().rho());
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // N_agg Equation
    {
        Info<< "N_agg Equation \n" << endl;
        tmp<fvScalarMatrix> N_aggEqn
        (
            fvm::ddt(rho, N_agg_)
            + fvm::div(phi, N_agg_)
            - fvm::laplacian(D*rho, N_agg_)
        ==
            - fvm::Sp(rho * S_coag_N_agg_ /N_agg_, N_agg_)
            + rho * S_inc_N_
        );

        N_aggEqn.ref().relax();
        fvOptions.constrain(N_aggEqn.ref());
        solve(N_aggEqn);
        fvOptions.correct(N_agg_);
    }

    // N_pri Equation
    {
        Info<< "N_pri Equation \n" << endl;
        tmp<fvScalarMatrix> N_priEqn
        (
            fvm::ddt(rho, N_pri_)
            + fvm::div(phi, N_pri_)
            - fvm::laplacian(D*rho, N_pri_)
        ==
            rho * S_inc_N_
        );

        N_priEqn.ref().relax();
        fvOptions.constrain(N_priEqn.ref());
        solve(N_priEqn);
        fvOptions.correct(N_pri_);
    }

    // C_tot Equation
    {
        Info<< "C_tot Equation \n" << endl;
        tmp<fvScalarMatrix> C_totEqn
        (
            fvm::ddt(rho, C_tot_)
            + fvm::div(phi, C_tot_)
            - fvm::laplacian(D*rho, C_tot_)
        ==
            rho * (S_inc_C_tot_ + S_grow_C_tot_ + S_ox_C_tot_)
        );

        C_totEqn.ref().relax();
        fvOptions.constrain(C_totEqn.ref());
        solve(C_totEqn);
        fvOptions.correct(C_tot_);
    }

    // H_tot Equation
    {
        Info<< "H_tot Equation \n" << endl;
        tmp<fvScalarMatrix> H_totEqn
        (
            fvm::ddt(rho, H_tot_)
            + fvm::div(phi, H_tot_)
            - fvm::laplacian(D*rho, H_tot_)
        ==
            rho * (S_inc_H_tot_ + S_grow_H_tot_)
        );

        H_totEqn.ref().relax();
        fvOptions.constrain(H_totEqn.ref());
        solve(H_totEqn);
        fvOptions.correct(H_tot_);
    }

}


// ************************************************************************* //
