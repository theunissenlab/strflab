function [modelParams, g] = glogGrad(modelParams, datIdx)
    
    global globDat;
    
    B = modelParams.B;
    M = modelParams.M;
    
    [modelParams, pr] = glogFwd(modelParams, datIdx);    
    pr = rv(pr);
    rd = rv(globDat.resp(datIdx)) - pr;
    
    hbm = exp(-B*(pr-M));    
    gbm = 1 + hbm;
    
    gb = (gbm.^-2) .* hbm .* (M - pr);
    gm = (gbm.^-2) .* hbm * B;
    
    bgrad = 2*rd .* gb;
    mgrad = 2*rd .* gm;
    
    g = [sum(mgrad) sum(bgrad)];