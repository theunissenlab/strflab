function [modelParams, g] = glogGradFD(modelParams, datIdx)

    [modelParams, w] = glogPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-8;
    
    
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = glogUnpak(modelParams, wpetFwd);
        modelParamsTempBack = glogUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = glogErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = glogErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    
    
    %{
    [modelParams, err] = glogErr(modelParams, datIdx);
    
    for k = 1:length(w)
        wpet = w;
        wpet(k) = w(k) + deps;
        
        modelParamsTemp = glogUnpak(modelParams, wpet);
        
        [modelParamsTemp, errFwd] = glogErr(modelParamsTemp, datIdx);
        
        g(k) = (errFwd - err) / deps;
    end
    %}