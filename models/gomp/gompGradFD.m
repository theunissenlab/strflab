function [modelParams, g] = gompGradFD(modelParams, datIdx)

    [modelParams, w] = gompPak(modelParams);
    g = zeros(1, length(w));
    
    deps = 1e-6;
    
    %{
    for k = 1:length(w)
        wpetBack = w;
        wpetFwd = w;
        wpetBack(k) = w(k) - deps;
        wpetFwd(k) = w(k) + deps;
        
        modelParamsTempFwd = gompUnpak(modelParams, wpetFwd);
        modelParamsTempBack = gompUnpak(modelParams, wpetBack);
        
        [modelParamsTempFwd, errFwd] = gompErr(modelParamsTempFwd, datIdx);
        [modelParamsTempBack, errBack] = gompErr(modelParamsTempBack, datIdx);
        
        g(k) = (errFwd - errBack) / (2*deps);
    end
    %}
    
    
    [modelParams, err] = gompErr(modelParams, datIdx);
    
    for k = 1:length(w)
        wpet = w;
        wpet(k) = w(k) + deps;
        
        modelParamsTemp = gompUnpak(modelParams, wpet);
        
        [modelParamsTemp, errFwd] = gompErr(modelParamsTemp, datIdx);
        
        g(k) = (errFwd - err) / deps;
    end
    