function modelParams = leglmInit(nIn, delays, family, dispersion)
    
    modelParams.type = 'leglm';
    
    modelParams.family = family;
    modelParams.nIn = nIn;
    modelParams.delays = delays;    
    modelParams.dispersion = dispersion;
    modelParams.w1 = zeros(nIn, length(delays));
    modelParams.b1 = 0;
    modelParams.m = 1;
    modelParams.nWts = nIn*length(delays) + 2;
    