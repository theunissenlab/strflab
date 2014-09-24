function [modelParams, resp] = gompFwd(modelParams, datIdx)

    global globDat;
    
    a = modelParams.a;
    b = modelParams.b;
    c = modelParams.c;
    
    resp = a*exp(b*exp(c*globDat.stim));
    resp = resp(datIdx);
    
    