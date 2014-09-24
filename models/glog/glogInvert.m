function [modelParams, invResp] = glogInvert(modelParams, resp)
    
    B = modelParams.B;
    M = modelParams.M;
    
    invResp = invglogistic(resp, B, M);
    