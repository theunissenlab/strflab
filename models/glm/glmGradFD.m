function [modelParams, g] = glmGradFD(modelParams, datIdx)

    [modelParams, w] = glmPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-8;
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = glmUnpak(modelParams, wpetFwd);
        modelParamsTempBack = glmUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = glmErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = glmErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    