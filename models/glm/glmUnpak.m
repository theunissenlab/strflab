function modelParams = glmUnpak(modelParams, w)

    nIn = modelParams.nIn;
    
    modelParams.w1 = reshape(w(1:nIn*length(modelParams.delays)), nIn, length(modelParams.delays));
    modelParams.b1 = w(end);
    
