function [modelParams, g] = gompGrad(modelParams, datIdx)
    
    global globDat;
    
    a = modelParams.a;
    b = modelParams.b;
    c = modelParams.c;
    
    [modelParams, pr] = gompFwd(modelParams, datIdx);    
    pr = rv(pr);
    rd = rv(globDat.resp(datIdx)) - pr;
    
    gc = b*exp(c*pr);
    bgrad = -a*exp(c*pr + gc);
    cgrad = -a*b*pr .* exp(c*pr + gc);
    
    bgrad = 2*rd .* bgrad;
    cgrad = 2*rd .* cgrad;
    
    g = [sum(bgrad) sum(cgrad)];
    
    