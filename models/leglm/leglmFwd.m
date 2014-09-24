function [modelParams, resp, linResp] = leglmFwd(modelParams, datIdx)

    global globDat;

    %% compute linear response
    linResp = conv_strf(globDat.stim(datIdx, :), modelParams.delays, modelParams.w1, globDat.groupIdx(datIdx)); %linear predictor
    linResp = linResp + modelParams.b1;
    
    resp = log(1 + exp(linResp)) .^ modelParams.m;    