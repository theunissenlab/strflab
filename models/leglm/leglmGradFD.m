function [modelParams, g] = leglmGradFD(modelParams, datIdx)

    [modelParams, w] = leglmPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-8;
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = leglmUnpak(modelParams, wpetFwd);
        modelParamsTempBack = leglmUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = leglmErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = leglmErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    