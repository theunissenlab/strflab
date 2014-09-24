function modelParams = npglmInit(nIn, delays, family, dispersion, basisFunc, centers)
    
    modelParams.type = 'npglm';
    
    modelParams.family = family;
    modelParams.nIn = nIn;
    modelParams.delays = delays;    
    modelParams.dispersion = dispersion;
    modelParams.w1 = zeros(nIn, length(delays));
    modelParams.b1 = 0;
    modelParams.nWts = nIn*length(delays) + 1 + length(centers);
    
    modelParams.basisFunc = basisFunc;
    modelParams.centers = centers;
    modelParams.w2 = zeros(1, length(centers));
    