function [modelParams, resp, respMat, dresp, R] = rbfFwd(modelParams, datIdx)

    global globDat;

    N = length(datIdx);
    M = modelParams.nfuncs;
    
    R = zeros(N, M);
    for k = 1:M        
        c = modelParams.centers(k, :);
        rtmp = bsxfun(@minus, globDat.stim(datIdx, :), c);        
        R(:, k) = sqrt(sum(rtmp, 2).^2);
    end
    
    [respMat, drespMat] = modelParams.func(R);
    dresp = drespMat * modelParams.w1';
    resp = respMat * modelParams.w1';
    