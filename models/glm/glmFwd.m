function [modelParams, resp, linResp] = glmFwd(modelParams, datIdx)

    global globDat;

    linResp = conv_strf(globDat.stim(datIdx, :), modelParams.delays, modelParams.w1, globDat.groupIdx(datIdx)); %linear predictor
    linResp = linResp + modelParams.b1;
    resp = modelParams.outputNL(linResp);
    