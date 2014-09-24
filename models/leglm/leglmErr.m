function [modelParams, err] = leglmErr(modelParams, datIdx)

    global globDat;

    [modelParams, u] = leglmFwd(modelParams, datIdx); % inverse link, mean of response
        
    y = globDat.resp(datIdx); %response variable
    a = modelParams.dispersion;     
    x = modelParams.family.canFunc(u); % canonical parameter
    b = modelParams.family.cumFunc(x); % cumulant function
    
    %{
    xnan = sum(isnan(x))
    bnan = sum(isnan(b))
    xinf = sum(isinf(x))
    binf = sum(isinf(b))
    %}
    
    err = -sum(y.*x - b) / a; %negative log likelihood
    