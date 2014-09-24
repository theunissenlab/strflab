function [modelParams, err, R] = rbfErr(modelParams, datIdx)

    global globDat;

    [modelParams, mresp, respMat, dresp, R] = rbfFwd(modelParams, datIdx);    
    d = globDat.resp(datIdx) - mresp;
    err = d' * d;    
    