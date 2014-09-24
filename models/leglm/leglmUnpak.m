function modelParams = leglmUnpak(modelParams, w)

    nIn = modelParams.nIn;
    
    w1len = nIn*length(modelParams.delays);
    modelParams.w1 = reshape(w(1:w1len), nIn, length(modelParams.delays));
    modelParams.b1 = w(w1len+1);
    modelParams.m = w(w1len+2);
    