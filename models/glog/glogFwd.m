function [modelParams, resp] = glogFwd(modelParams, datIdx)

    global globDat;
    
    M = modelParams.M;
    B = modelParams.B;
    
    resp = glogistic(globDat.stim, B, M);
    resp = resp(datIdx);
    