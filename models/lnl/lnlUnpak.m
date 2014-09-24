function modelParams = lnlUnpak(modelParams, w)

    nlin = modelParams.linModel.nWts;
    
    wlin = w(1:nlin);
    wnl = w(nlin+1:end);
    
    modelParams.linModel = strfUnpak(modelParams.linModel, wlin);
    modelParams.nlModel = strfUnpak(modelParams.nlModel, wnl);
    