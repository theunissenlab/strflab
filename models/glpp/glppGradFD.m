function [modelParams, g] = glppGradFD(modelParams, datIdx)

    [modelParams, w] = glppPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-10;
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = glppUnpak(modelParams, wpetFwd);
        modelParamsTempBack = glppUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = glppErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = glppErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    
    