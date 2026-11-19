// Called once before each hybrid PTC update with the preceding steady residual.
// Store the detector and trial state in the case's checkpoint dictionary.
const dictionary feedback = controls.subOrEmptyDict("pseudoCoResidualFeedback");
if (feedback.lookupOrDefault<bool>("enabled", false))
{
    const scalar targetDecay = feedback.lookupOrDefault<scalar>
    (
        "targetDecay", log(scalar(10))*8/8000
    );
    const scalar upThreshold = feedback.lookupOrDefault<scalar>("upThreshold", 0.15);
    const scalar downThreshold = feedback.lookupOrDefault<scalar>
    (
        "downThreshold", 100*targetDecay
    );
    const scalar increase = feedback.lookupOrDefault<scalar>("increase", 1.02);
    const scalar decrease = feedback.lookupOrDefault<scalar>("decrease", 0.8);
    const scalar minFactor = feedback.lookupOrDefault<scalar>("minFactor", 1e-4);
    const scalar maxFactor = feedback.lookupOrDefault<scalar>("maxFactor", 100);
    const scalar referenceWeight = feedback.lookupOrDefault<scalar>("referenceWeight", 0.1);
    const scalar trialTolerance = feedback.lookupOrDefault<scalar>("trialTolerance", 0.15);
    const label trialMaxSteps = feedback.lookupOrDefault<label>("trialMaxSteps", 40);
    for (const scalar value : {targetDecay, upThreshold, downThreshold, increase,
        decrease, minFactor, maxFactor, referenceWeight, trialTolerance})
    {
        if (!(std::isfinite(value) && value > 0))
            FatalErrorInFunction << "Residual CFL feedback settings must be finite and positive"
                << exit(FatalError);
    }
    if (increase <= 1 || decrease >= 1 || minFactor > 1 || maxFactor < 1
        || referenceWeight > 1 || trialMaxSteps < 1)
        FatalErrorInFunction << "Invalid residual CFL feedback factors, weight or trial length"
            << exit(FatalError);

    IOdictionary& history = mesh_.lookupObjectRef<IOdictionary>("hisaCFLFeedbackState");
    const label index = mesh_.time().timeIndex();
    const label lastIndex = history.lookupOrDefault<label>("lastIndex", -2);
    const scalar logResidual = log(relativeBounded);
    const scalar currentCFL = pseudoCoNum_().value();
    const scalar baseCFL = safety*initial/denominator;
    scalar factor = history.lookupOrDefault<scalar>("factor", 1);
    scalar upIntegral = history.lookupOrDefault<scalar>("upIntegral", 0);
    scalar downIntegral = history.lookupOrDefault<scalar>("downIntegral", 0);
    scalar meanDecay = history.lookupOrDefault<scalar>("meanDecay", 0);
    scalar previousFactor = history.lookupOrDefault<scalar>("previousFactor", 1);
    scalar trialStart = history.lookupOrDefault<scalar>("trialStartLogResidual", logResidual);
    scalar trialReference = history.lookupOrDefault<scalar>("trialReferenceDecay", 0);
    label trialSteps = history.lookupOrDefault<label>("trialSteps", -1);
    bool waitForTrend = history.lookupOrDefault<bool>("waitForTrend", false);
    word event("none");

    if (lastIndex == index-1)
    {
        const scalar decay = history.get<scalar>("logResidual") - logResidual;
        const scalar error = decay-targetDecay;
        meanDecay += referenceWeight*(decay-meanDecay);
        upIntegral = max(scalar(0), upIntegral+error);
        downIntegral = max(scalar(0), downIntegral-error);

        if (trialSteps >= 0)
        {
            ++trialSteps;
            const scalar improvement = trialStart-logResidual-trialSteps*trialReference;
            history.set("trialImprovement", improvement);
            if (improvement >= trialTolerance)
            {
                trialSteps = -1;
                event = "acceptTrial";
            }
            else if (improvement <= -trialTolerance || trialSteps >= trialMaxSteps)
            {
                factor = previousFactor;
                trialSteps = -1;
                waitForTrend = true;
                upIntegral = downIntegral = 0;
                event = "rejectTrial";
            }
        }
        else
        {
            if (mag(meanDecay-trialReference) >= trialTolerance/trialMaxSteps)
                waitForTrend = false;

            const scalar trialFactor = max(minFactor, factor*decrease);
            if (!waitForTrend && downIntegral >= downThreshold
                && currentCFL > pseudoCoNumMin_ && trialFactor < factor
                && baseCFL*trialFactor < currentCFL)
            {
                previousFactor = factor;
                factor = trialFactor;
                trialStart = logResidual;
                trialReference = meanDecay;
                trialSteps = 0;
                downIntegral -= downThreshold;
                upIntegral = 0;
                history.set("trialImprovement", scalar(0));
                event = "startTrial";
            }
            else if (upIntegral >= upThreshold && factor < maxFactor
                && currentCFL < pseudoCoNumMax_)
            {
                factor = min(maxFactor, factor*increase);
                upIntegral -= upThreshold;
                event = "increase";
            }
        }
    }
    else if (lastIndex != index)
    {
        // A gap in PTC (e.g. an NK excursion) invalidates a trial comparison.
        if (trialSteps >= 0)
            factor = previousFactor;
        trialSteps = -1;
        waitForTrend = false;
        upIntegral = downIntegral = meanDecay = 0;
        event = "initial";
    }

    factor = max(minFactor, min(maxFactor, factor));
    const scalar targetCFL = baseCFL*factor;
    ratio = targetCFL/currentCFL;
    history.set("lastIndex", index);
    history.set("logResidual", logResidual);
    history.set("factor", factor);
    history.set("upIntegral", upIntegral);
    history.set("downIntegral", downIntegral);
    history.set("meanDecay", meanDecay);
    history.set("previousFactor", previousFactor);
    history.set("trialStartLogResidual", trialStart);
    history.set("trialReferenceDecay", trialReference);
    history.set("trialSteps", trialSteps);
    history.set("waitForTrend", waitForTrend);
    if (event != "none")
    {
        history.set("lastEvent", event);
        history.set("lastEventIndex", index);
    }
    history.set("targetCFL", targetCFL);
}
