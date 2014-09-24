function [modelParams, g] = rbfGradFD(modelParams, datIdx)

    [modelParams, w] = rbfPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-8;
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = rbfUnpak(modelParams, wpetFwd);
        modelParamsTempBack = rbfUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = rbfErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = rbfErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    