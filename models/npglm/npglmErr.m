function [modelParams, err] = npglmErr(modelParams, datIdx)

    global globDat;

    [modelParams, u] = npglmFwd(modelParams, datIdx); % inverse link, mean of response
        
    y = globDat.resp(datIdx); %response variable
    a = modelParams.dispersion;     
    x = modelParams.family.canFunc(u); % canonical parameter
    b = modelParams.family.cumFunc(x); % cumulant function

    err = -sum(y.*x - b) / a; %negative log likelihood
    