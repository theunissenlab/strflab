function [modelParams, g] = glppGradLowMem(modelParams, datIdx)

    global globDat;

    if modelParams.useFDGrad
        %fprintf('Using finite difference approximation to gradient.\n');
        [modelParams, gfd] = glppGradFD(modelParams, datIdx);
	gfd
        g = gfd;
        return;
    end
    
    M = modelParams.nIn;
    P = length(modelParams.delays);
    nStimWts = M*P + ~modelParams.ignoreBias;
    nSpikeWts = length(modelParams.spikeResponseWeights);
    f = @(x) modelParams.nlFunc(x);
    df = @(x) modelParams.nlFuncDerivative(x);
    swts = modelParams.spikeResponseWeights;
    numTrials = length(globDat.resp);
    dt = 1 / modelParams.sampleRate;
    minLogTol = 1e-10;
    
    %% gradient is composed of 4 parts, g1{Stim,Spike} is the gradient
    %% of the first term in the likelihood (L1) with response to stim and
    %% spike weights. g2{Stim,Spike} is the gradient of L2
    g1Stim = zeros(1, nStimWts);
    g1Spike = zeros(1, nSpikeWts);
    
    g2Stim = zeros(1, nStimWts);
    g2Spike = zeros(1, nSpikeWts);
    
    stimAtT = zeros(P, modelParams.nIn);
    
    %% compute stimulus current
    stim = globDat.stim(datIdx);
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
    bias = 0;
    if ~modelParams.ignoreBias
        bias = modelParams.b1;
    end
    stimCurrent = conv_strf(stim, modelParams.delays, modelParams.w1, groupIndex) + bias;
    
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
                if t <= length(spikeCurrent)
                    eindx = min(length(spikeCurrent), t+nSpikeWts);
                    rng = (t+1):eindx;

                    spikeCurrent(rng) = spikeCurrent(rng) + swts(1:length(rng));
                end

            end
        end

        totalCurrent = stimCurrent + spikeCurrent;
        rate = f(totalCurrent);
        drate = df(totalCurrent);

        %% compute derivatives for first and second component of error term
        for t = 1:size(stim, 1)

            isSpike = ~isempty(find(spikeIndicies == t));

            r = rate(t);
            dr = drate(t);
            if r <= 0
                error('Rate function for lnp model must be non-negative.\n');
            end

            if ~modelParams.fixedStrf
                sindx = t - P + 1;
                sindx = max(sindx, 1);
                eindx = max(P - t + 1, 1);

                %compute stimulus driven part of gradient
                stimAtT(:, :) = 0;
                stimAtT(eindx:P, :) = globDat.stim(sindx:t, :);

                istim = flipud(stimAtT);
                istim = istim';
                istim = [istim(:)'];
                if ~modelParams.ignoreBias
                    istim = [istim 1];
                end

                if isSpike
                    gval = (dr/r)*istim;
                    g1Stim = g1Stim + gval;
                end

                g2Stim = g2Stim + dr*istim;
            end

            %compute spike response part of gradient
            if ~isempty(swts)
                gspike = zeros(1, length(swts));
                prevSpikeIntervals = t - spikeIndicies(spikeIndicies < t);
                prevSpikeIntervals = prevSpikeIntervals(prevSpikeIntervals <= length(swts));
                for m = 1:length(swts) 
                    mindx = find(prevSpikeIntervals == m);
                    if ~isempty(mindx)
                        gspike(m) = length(mindx);
                    end
                end

                if isSpike
                    g1Spike = g1Spike + (gspike * (dr/r));
                end
                g2Spike = g2Spike + gspike*dr;
            end

        end

    end
        
    g1 = [g1Stim g1Spike];
    g2 = [g2Stim g2Spike]*dt; %g2 is an integral, make sure it's multiplied by dt
   
    %% compute derivative for regularizer
    gr = zeros(1, nStimWts + nSpikeWts);
    if modelParams.regularize
        Rwt = modelParams.regularizeWeight;
        gr = Rwt*2*[modelParams.w1(:)' modelParams.b1 modelParams.spikeResponseWeights];
    end

    if modelParams.debugL1
        fprintf('Leaving out L2..\n');
        g2(:) = 0;
        gr = 0;
    end
    if modelParams.debugL2
        fprintf('Leaving out L1..\n');
        g1(:) = 0;
        gr = 0;
    end
    g = -(g1 - g2 - gr);

    if modelParams.checkGrad
        [modelParams, gfd] = glppGradFD(modelParams, datIdx);
        g
        gfd	
        gdiff = norm(gfd - g);
        ndiff = gdiff / (0.5*(norm(gfd)+norm(g)));
        fprintf('Gradient check, diff=%f, ndiff=%f\n', gdiff, ndiff);
        %g
        %gfd
    end
    