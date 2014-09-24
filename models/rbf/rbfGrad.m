function [modelParams, g] = rbfGrad(modelParams, datIdx)

    global globDat;    

    [modelParams, mresp, respMat, dresp, R] = rbfFwd(modelParams, datIdx);    
    d = globDat.resp(datIdx) - mresp;    
    pR = modelParams.func(R');    
    g = -2*(pR * d)';    
    