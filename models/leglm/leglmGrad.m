function [modelParams, g] = leglmGrad(modelParams, datIdx)

    global globDat;    

    stim = globDat.stim(datIdx, :);
    groupIndex = globDat.groupIdx(datIdx);
    
    [modelParams, u, linResp] = leglmFwd(modelParams, datIdx);
    
    m = modelParams.m;
    
    y = globDat.resp(datIdx); %response variable    
    a = modelParams.dispersion;
    iv = modelParams.family.meanVarInv(u); % inverse of variance of mean

    d = y - u;
    div = d .* iv;
    elr = exp(linResp);
    lbase = log(1+elr);
    dle = m*(lbase.^(m-1));
    dg = dle .* (elr ./ (1 + elr));
    r = div .* dg;
    
    stimMat = zeros(length(modelParams.delays), modelParams.nIn);
    gSum = zeros(length(modelParams.delays), modelParams.nIn);    
    
    grps = unique(groupIndex);
    for k = 1:length(grps)
        gindx = find(groupIndex == grps(k));        
        gstart = min(gindx);
        for t = gindx
            stimMat(:, :) = 0;
            soff = t - length(modelParams.delays) + 1;
            s = max(gstart, soff);
            me = length(modelParams.delays) - (t-s);
            
            stimMat(me:end, :) = stim(s:t, :);
            gSum = gSum + (stimMat*r(t));
            
        end
    end
    
    gSum = fliplr(gSum');
    gw = gSum(:)';
    gb = sum(r);    
    gm = div * (u .* log(lbase))';
            
    g = -[gw gb gm] / a;
    