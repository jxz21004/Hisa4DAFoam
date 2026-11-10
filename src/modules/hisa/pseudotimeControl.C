/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2015 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2015 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2022 Oliver Oxtoby
    Copyright (C) 2011-2013 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pseudotimeControl.H"
#include "Switch.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudotimeControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

#if FOUNDATION >= 6 or OPENFOAM >= 1906
bool Foam::pseudotimeControl::read()
#else
void Foam::pseudotimeControl::read()
#endif
{
#if FOUNDATION >= 6
    singleRegionSolutionControl::read();
#else
    solutionControl::read(false);
#endif

    // Read solution controls
    const dictionary& pseudoDict = dict();
    nCorrOuter_ = pseudoDict.lookupOrDefault<label>("nPseudoCorr", 20);
    nCorrOuterMin_ = pseudoDict.lookupOrDefault<label>("nPseudoCorrMin", 1);
    if (pseudoDict.found("pseudoTol"))
    {
        pseudoDict.lookup("pseudoTol") >> residualTols_;
    }
    if (pseudoDict.found("pseudoTolRel"))
    {
        pseudoDict.lookup("pseudoTolRel") >> residualTolsRel_;
    }
    turbOnFinalIterOnly_ =
        pseudoDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", false);

    if (state_.found("initResiduals"))
    {
        state_.lookup("initResiduals") >> initResiduals_;
    }
    else
    {
        initResiduals_ = residualIO(residualTols_, 0.0);
    }
    firstResiduals_.clear();
    if (state_.found("firstResiduals"))
    {
        firstResiduals_.reset(new residualIO(residualTols_, 0.0));
        state_.lookup("firstResiduals") >> firstResiduals_();
    }
#if FOUNDATION >= 6 or OPENFOAM >= 1906
    return true;
#endif
}


bool Foam::pseudotimeControl::criteriaSatisfied()
{
    // no checks on first iteration after restart or first iteration of
    // timestep - nothing has been calculated yet
    if (firstIteration_ || corr_ == 1 || finalIter())
    {
        firstIteration_ = false;
        return false;
    }

    bool storeIni = this->storeInitialResiduals();

    if (storeIni)
    {
        initResiduals_ = residuals_;
        state_.add("initResiduals", initResiduals_, true); // Store for restarts
        if (!firstResiduals_.valid())
        {
           firstResiduals_.reset(new residualIO(residuals_));
           state_.add("firstResiduals", firstResiduals_(), true); // Store for restarts
        }
    }

    bool absCheck = (max(residuals_-residualTols_) < 0.0);

    bool relCheck = false;

    if (!storeIni)
    {
        relCheck = (max(residuals_/(initResiduals_+ROOTVSMALL) - residualTolsRel_) < 0.0);
    }

    return (corr_ >= nCorrOuterMin_) && (absCheck || relCheck);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bool Foam::pseudotimeControl::checkFixedFlowResidual(const scalarField& rms)
{
    if (!dict().lookupOrDefault<Switch>("fixedInitialResidual", false)) return false;
    if (!steadyState_) FatalErrorInFunction << "fixedInitialResidual requires steady flow" << exit(FatalError);
    scalarField initial;
    const bool fresh = !state_.found("fixedFlowInitialRMS");
    if (fresh)
    {
        if (mesh_.time().startTime().value() != 0 && !dict().lookupOrDefault<Switch>("allowNewResidualBaseline", false))
            FatalErrorInFunction << "Restart has no fixedFlowInitialRMS; supply original baseline or explicitly allowNewResidualBaseline true" << exit(FatalError);
        initial=rms;
        state_.add("fixedFlowInitialRMS", initial);
    }
    else state_.lookup("fixedFlowInitialRMS") >> initial;
    if (initial.size()!=5 || rms.size()!=5) FatalErrorInFunction << "Expected five flow residuals" << exit(FatalError);
    scalarField absolute(5,0.0);
    if (dict().found("fixedResidualAbsTol")) dict().lookup("fixedResidualAbsTol") >> absolute;
    const scalar tol=dict().lookupOrDefault<scalar>("fixedResidualRelTol",1e-8);
    if (absolute.size()!=5 || !std::isfinite(tol) || tol<=0 || tol>=1) FatalErrorInFunction << "Invalid fixed residual tolerances" << exit(FatalError);
    scalarField relative(5,0.0);
    bool pass=true;
    forAll(rms,i)
    {
        if (!std::isfinite(rms[i]) || !std::isfinite(initial[i]) || rms[i]<0 || initial[i]<0 || !std::isfinite(absolute[i]) || absolute[i]<0)
            FatalErrorInFunction << "Nonfinite or invalid steady residual/baseline/tolerance" << exit(FatalError);
        relative[i]=initial[i]>0 ? rms[i]/initial[i] : (rms[i]==0 ? 0 : GREAT);
        pass = pass && (rms[i]<=absolute[i] || (initial[i]>0 && relative[i]<=tol));
    }
    Info<< "FixedFlowResidualRMS: " << mesh_.time().value();
    forAll(rms,i) Info<< " " << rms[i];
    Info<< nl << "FixedFlowResidualRelative: " << mesh_.time().value();
    forAll(relative,i) Info<< " " << relative[i];
    Info<< nl << "FixedFlowResidualInitial:";
    forAll(initial,i) Info<< " " << initial[i];
    Info<< endl;
    fixedFlowConverged_=!fresh && corr_>=nCorrOuterMin_ && pass;
    if(fixedFlowConverged_) Info<< "Fixed flow residual criterion satisfied (SA not included)" << endl;
    return fixedFlowConverged_;
}

Foam::pseudotimeControl::pseudotimeControl(fvMesh& mesh, const bool steadyState, const label nScalars, const label nVectors, const residualIO& defaultTol, const residualIO& defaultTolRel)
:
#if FOUNDATION >= 6
    singleRegionSolutionControl(mesh, "pseudoTime"),
    corr_(0),
#else
    solutionControl(mesh, "pseudoTime"),
#endif
    steadyState_(steadyState),
    nCorrOuter_(0),
    turbOnFinalIterOnly_(false),
    converged_(false),
    firstIteration_(true),
    residualTols_(defaultTol),
    residualTolsRel_(defaultTolRel),
    residuals_(defaultTol, 0.0),
    initResiduals_(defaultTol, 0.0),
    state_
    (
        IOobject
        (
            "pseudotimeState",
            #if FOUNDATION >= 12
            mesh.time().name(),
            #else
            mesh.time().timeName(),
            #endif
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    read();

    Info<< algorithmName_ << ": ";
    if (!steadyState)
    {
        Info<< "max iterations = " << nCorrOuter_ << ", ";
    }
    Info<< "tolerance = " << residualTols_ << ", "
        << "relTol = " << residualTolsRel_
        << nl;
    Info<< endl;

    printInterval_ = mesh.time().controlDict().lookupOrDefault<label>("HisaPrintInterval", 1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudotimeControl::~pseudotimeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pseudotimeControl::loop()
{
    read();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    // If steady state, ignore nCorrOuter because endTime is used
    // to terminate main loop
    if (!steadyState_ && corr_ == nCorrOuter_ + 1)
    {
        if (nCorrOuter_ != 1)
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrOuter_ << " iterations" << endl;
        }

        corr_ = 0;
        if (!steadyState_)
        {
                #if FOUNDATION >= 12
                finalIter_ = false;
                #elif OPENFOAM >= 2312
            mesh_.data().setFinalIteration(false);
            #else
            mesh_.data::remove("finalIteration");
            #endif
        }
        return false;
    }

    if (steadyState_ && dict().lookupOrDefault<Switch>("fixedInitialResidual", false))
    {
        if (fixedFlowConverged_) return false;
        // Preserve legacy initial-residual bookkeeping used by CFL adaptation,
        // but never use its dynamic normalization as the stopping criterion.
        (void)criteriaSatisfied();
        Info<< algorithmName_ << ": iteration " << corr_ << endl;
        return true;
    }
    bool completed = false;
    if (converged_) converged_ = (max(residuals_-residualTols_) < 0.0)
        || (max(residuals_/(initResiduals_+ROOTVSMALL)-residualTolsRel_) < 0.0);
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;

            if (!steadyState_)
            {
                #if FOUNDATION >= 12
                finalIter_ = false;
                #elif OPENFOAM >= 2312
                mesh_.data().setFinalIteration(false);
                #else
                mesh_.data::remove("finalIteration");
                #endif
            }
            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();

            if (!steadyState_)
            {
                #if FOUNDATION >= 12
                finalIter_ = true;
                #elif OPENFOAM >= 2312
                mesh_.data().setFinalIteration(true);
                #else
                mesh_.data::add("finalIteration", true);
                #endif
            }
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            if (!steadyState_)
            {
                #if FOUNDATION >= 12
                finalIter_ = true;
                #elif OPENFOAM >= 2312
                mesh_.data().setFinalIteration(true);
                #else
                mesh_.data::add("finalIteration", true);
                #endif
            }
        }

        if (steadyState_ || corr_ <= nCorrOuter_)
        {
            if (nCorrOuter_ != 1)
            {
                if (mesh_.time().timeIndex() % printInterval_ == 0 || this->mesh_.time().timeIndex() == 1)
                {
                    Info<< algorithmName_ << ": iteration " << corr_ << endl;
                }

                storePrevIterFields();
            }

            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
