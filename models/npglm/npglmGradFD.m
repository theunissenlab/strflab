function [modelParams, g] = npglmGradFD(modelParams, datIdx)

    [modelParams, w] = npglmPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-8;
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = npglmUnpak(modelParams, wpetFwd);
        modelParamsTempBack = npglmUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = npglmErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = npglmErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    