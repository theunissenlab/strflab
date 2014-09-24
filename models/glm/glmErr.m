function [modelParams, err] = glmErr(modelParams, datIdx)

    global globDat;

    linResp = conv_strf(globDat.stim(datIdx, :), modelParams.delays, modelParams.w1, globDat.groupIdx(datIdx)); %linear predictor
    linResp = linResp + modelParams.b1;
    y = globDat.resp(datIdx); %response variable

    u = modelParams.outputNL(linResp); % inverse link, mean of response
    a = modelParams.dispersion;     
    x = modelParams.family.canFunc(u); % canonical parameter
    b = modelParams.family.cumFunc(x); % cumulant function

    err = -sum(y.*x - b) / a; %negative log likelihood
    
    %% regularization
    if modelParams.qwt > 0
        w = modelParams.w1(:);
        r = (w.^2) .^ (1/qval);
        err = err + qwt*sum(r);
    end
