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
    coagulation_enabled_(true)
    // PAH_names_(
    //     sootProps_.get<wordList>("PAHs")
    // ),
    // PAH_n_C_(PAH_names_.size()),
    // PAH_n_H_(PAH_names_.size()),
    // PAH_indicies_(PAH_names_.size()),
    // dimer_names_(PAH_names_.size() * (PAH_names_.size() +1 ) / 2),
    // dimer_n_C_(dimer_names_.size()),
    // dimer_n_H_(dimer_names_.size()),
    // dimer_PAH_1_index_(dimer_names_.size()),
    // dimer_PAH_2_index_(dimer_names_.size()),
    // dimer_PAH_1_id_(dimer_names_.size()),
    // dimer_PAH_2_id_(dimer_names_.size())
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

// Return Molecular Weight in kg/mol
// template<class ReactionThermo>
// dimensionedScalar Foam::combustionModels::laminarSoot<ReactionThermo>::W(label index) const
// {
//     basicSpecieMixture& composition = this->thermo().composition();
//     return composition.W(index) / 1000 * dimensionedScalar("onekgPerMol", dimensionSet(1,0,0,0,-1,0,0), scalar(1));
// }

// Return concentration in mol/m3
// template<class ReactionThermo>
// volScalarField Foam::combustionModels::laminarSoot<ReactionThermo>::C(label index) const
// {
//     basicSpecieMixture& composition = this->thermo().composition();
//     return composition.Y()[index] * this->thermo().rho() / W(index);
// }

// ************************************************************************* //
