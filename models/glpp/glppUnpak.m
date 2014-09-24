function modelParams = glppUnpak(modelParams, w)

    nIn = modelParams.nIn;
    strfEndIndex = nIn*length(modelParams.delays);
    modelParams.w1 = reshape(w(1:strfEndIndex), nIn, length(modelParams.delays));
    
    if ~modelParams.ignoreBias
        modelParams.b1 = w(strfEndIndex + 1);
        strfEndIndex = strfEndIndex + 1;
    end

    modelParams.spikeResponseWeights = w((strfEndIndex+1):end);
    