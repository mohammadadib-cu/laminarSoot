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
    coagulation_enabled_(true),
    PAH_names_(
        sootProps_.get<wordList>("PAHs")
        // sootProps_.lookup("PAHs")
    ),
    PAH_n_C_(PAH_names_.size()),
    PAH_n_H_(PAH_names_.size()),
    PAH_indicies_(PAH_names_.size())

{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }

    readPAHs();
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


// ************************************************************************* //
