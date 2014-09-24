function [modelParams, g] = npglmGrad(modelParams, datIdx)

    global globDat;    

    stim = globDat.stim(datIdx, :);
    groupIndex = globDat.groupIdx(datIdx);
    
    [modelParams, u, ~, Rdiff, uMat, duMat] = npglmFwd(modelParams, datIdx);
    
    y = globDat.resp(datIdx); %response variable    
    a = modelParams.dispersion;
    iv = modelParams.family.meanVarInv(u); % inverse of variance of mean

    d = y - u;
    div = d .* iv;
    
    Rsig = sign(Rdiff);
    
    nt = (Rsig .* duMat) * modelParams.w2';
    
    r = div .* nt;
    
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
    
    gnl = uMat' * div;
    
    g = -[gw gb gnl'] / a;
    