function [modelParams, g] = glmGrad(modelParams, datIdx)

    global globDat;    

    stim = globDat.stim(datIdx, :);
    groupIndex = globDat.groupIdx(datIdx);
    
    linResp = conv_strf(stim, modelParams.delays, modelParams.w1, globDat.groupIdx(datIdx)); %linear predictor
    linResp = linResp + modelParams.b1;
    y = globDat.resp(datIdx); %response variable
    
    [u, du] = modelParams.outputNL(linResp); % inverse link, mean of response
    a = modelParams.dispersion;
    iv = modelParams.family.meanVarInv(u); % inverse of variance of mean

    r = (y - u) .* iv .* du;
    
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
    gb = sum(r);
    
    g = -[gSum(:)' gb] / a;