function [modelParams, resp, R, Rdiff, respMat, drespMat] = npglmFwd(modelParams, datIdx)

    global globDat;

    %% compute linear response
    linResp = conv_strf(globDat.stim(datIdx, :), modelParams.delays, modelParams.w1, globDat.groupIdx(datIdx)); %linear predictor
    linResp = linResp + modelParams.b1;
    
    N = length(datIdx);
    M = length(modelParams.centers);
    
    Rdiff = zeros(N, M);    
    for k = 1:M        
        c = modelParams.centers(k, :); %should be 1D (not built to handle higher dimensions yet...)
        Rdiff(:, k) = linResp - c;
    end
    
    R = abs(Rdiff); %1D shortcut, should be sqrt(rdiff.^2)
    
    [respMat, drespMat] = modelParams.basisFunc(R);    
    resp = respMat * modelParams.w2';
    