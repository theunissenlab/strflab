function [modelParams, err] = glppErr(modelParams, datIdx)    
    global globDat;

    %% initialize some preliminary stuff
    srLen = length(modelParams.spikeResponseWeights);
    swts = modelParams.spikeResponseWeights;
    f = @(x) modelParams.nlFunc(x);
    numTrials = length(globDat.resp);
    dt = 1 / modelParams.sampleRate;
    
    L1 = 0;
    L2 = 0;
    
    %% compute stimulus current once, same across trials
    stim = globDat.stim;
    %rearrange group indicies to accomodate held-out trials
    groupIndex = globDat.groupIdx;
    if ~isempty(datIdx)
        stim = globDat.stim(datIdx, :);
        groupIndex = globDat.groupIdx(datIdx);
        uindx = unique(groupIndex);

        for k = 1:length(uindx)
            ui = uindx(k);
            groupIndex(groupIndex == ui) = k;
        end

    end
    %add bias
    bias = 0;
    if ~modelParams.ignoreBias
        bias = modelParams.b1;
    end
    %convolve filter with stimulus
    stimCurrent = conv_strf(stim, modelParams.delays, modelParams.w1, groupIndex) + bias;
   
    %% compute error terms for each trial
    for trialNum = 1:numTrials
        
        %get spike times in seconds
        spikeTimes = globDat.resp{trialNum};
        %convert spike times to stimulus indicies
        spikeIndicies = round(spikeTimes*modelParams.sampleRate) + 1;
        
        %% compute spike response current
        spikeCurrent = zeros(size(stimCurrent));
        if ~isempty(swts)
            for k = 1:length(spikeIndicies)

                t = spikeIndicies(k);
                eindx = min(length(spikeCurrent), t+srLen);
                rng = (t+1):eindx;

                spikeCurrent(rng) = spikeCurrent(rng) + swts(1:length(rng));

            end
        end

        totalCurrent = stimCurrent + spikeCurrent;
        rate = f(totalCurrent);
        
        %% compute L1 for this trial
        for k = 1:length(spikeIndicies)

            t = spikeIndicies(k);
            
            if t <= length(rate)                
                rval = rate(t);
                if rval <= 0                    
                    %rval = 1e-3;
                    error('Rate function in lnp model must be non-negative, rval=%f\n', rval);
                end
                L1 = L1 + log(rval);

                %check for problems
                if isinf(L1)
                    sval = totalCurrent(t);
                    strf = modelParams.w1
                    bias
                    swts
                    warning('Error=Inf, try perturbing initial guess and rerunning:\nL1=Inf, sval=%f, rval=%f\n', sval, rval);
                end
            end
        end

        %% compute L2 for this trial
        L2 = L2 + trapz(rate)*dt;
        
    end
    
    %fprintf('stimSum=%f, spikeSum=%f, rateSum=%f\n', sum(stimCurrent), sum(spikeCurrent), sum(rate));
    
    %% compute regularization
    R = 0;
    if modelParams.regularize
        lambda = modelParams.lambda;
        [modelParams, w] = glppPak(modelParams);
        r1 = (1 - lambda)*sum(abs(w));
        r2 = lambda*sum(w.^2);
        R = modelParams.normC*(r1 + r2);
    end

    if modelParams.debugL1
        fprintf('Leaving out L2..\n');
        L2 = 0;
        R = 0;
    end
    if modelParams.debugL2
        fprintf('Leaving out L1..\n');
        L1 = 0;
        R = 0;
    end
    
    %% error is negative regularized log likelihood
    likelihood = L1 - L2;
    %if isfield(modelParams, 'normC')
    %    likelihood = likelihood / modelParams.normC;
    %end
    err = -(likelihood - R);
    
    %strf = modelParams.w1
    %bias = modelParams.b1
    %swts = modelParams.spikeResponseWeights
    fprintf('L1=%f, L2=%f, R=%f, err=%f\n', L1, L2, R, err);
    