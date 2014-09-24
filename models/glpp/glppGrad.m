function [modelParams, g] = glppGrad(modelParams, datIdx)

    global globDat;

    if modelParams.useFDGrad
        %fprintf('Using finite difference approximation to gradient.\n');
        [modelParams, gfd] = glppGradFD(modelParams, datIdx);
        g = gfd;
        return;
    end
    
    if modelParams.useLowMemory
        %fprintf('Using low-memory gradient (slower for high dimensional stim).\n');
        [modelParams, gLM] = glppGradLowMem(modelParams, datIdx);
        g = gLM;
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
    
    %% gradient is composed of 4 parts, g1{Stim,Spike} is the gradient
    %% of the first term in the likelihood (L1) with response to stim and
    %% spike weights. g2{Stim,Spike} is the gradient of L2
    ngStimWts = nStimWts - ~modelParams.ignoreBias;
    g1Stim = zeros(1, ngStimWts);
    g1Spike = zeros(1, nSpikeWts);
    
    g2Stim = zeros(1, ngStimWts);
    g2Spike = zeros(1, nSpikeWts);
    
    stimAtT = zeros(P, modelParams.nIn);
    
    %% compute stimulus current
    stim = globDat.stim(datIdx);
    %rearrange group indicies so they go 1,2,3,... (no gaps from left out stimuli)
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
        
    stimLen = size(stim, 1);    
    if ~modelParams.fixedStrf
        %pre-allocate stimulus history matrix
        
        maxMem = modelParams.maxGradMem; %max memory for stimulus history matrix in MB
        memMult = 8 / (1024*1024); %MB per element
        maxElems = round(maxMem / (memMult*M));
        histChunkLen = min(maxElems, stimLen);
        
        shist = zeros(ngStimWts, histChunkLen);
        stimAtT = zeros(P, modelParams.nIn);
        %fprintf('stimulus history matrix takes up around %0.2fMB\n', maxMem);
    end
    
    g1Bias = 0;
    g2Bias = 0;
    
    %% go through each trial and compute it's contribution to the gradient
    for trialNum = 1:numTrials
    
        %get spike times in seconds
        spikeTimes = globDat.resp{trialNum};
        %convert spike times to stimulus indicies
        spikeIndicies = round(spikeTimes*modelParams.sampleRate) + 1;
        spindx = spikeIndicies;
        
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
        irate = rate.^-1;
        rhs = (irate(spindx) .* drate(spindx));

        %% form matricies of stimulus history vectors and compute stimulus
        %  driven gradient in chunks
        %tic;
        if ~modelParams.fixedStrf
            
            numChunks = ceil(stimLen / histChunkLen);
            for k = 1:numChunks
                
                shist(:, :) = 0;
                tstart = (k-1)*histChunkLen + 1;
                tend = min(tstart + histChunkLen, stimLen);
                trange = tstart:tend;
                
                for m = 1:length(trange)
                    t = trange(m);
                    %% get stimulus history at time t, change it to a vector,
                    %  flip it, and put it in the matrix
                    sindx = max(t-P+1, 1);
                    eindx = max(P - t + 1, 1);
                    stimAtT(:, :) = 0;
                    stimAtT(eindx:P, :) = globDat.stim(sindx:t, :);            
                    istim = flipud(stimAtT);
                    istim = istim';
                    istim = istim(:)';

                    %throw it in matrix
                    shist(:, m) = istim;

                end
                
                %% do matrix multiplications to get stimulus-driven part of gradient
                sspindx = spindx( (spindx >= min(trange)) & (spindx <= max(trange)) );
                srhs = (irate(sspindx) .* drate(sspindx));
                g1Add = (shist(:, sspindx) * srhs');
                g1Stim = g1Stim + g1Add';

                g2Add = (shist * drate(trange)');
                g2Stim = g2Stim + g2Add';
                
            end
                
        end
        %tstim = toc;
        
        %% form elements for sparse matrix of indicator vectors for gradient w/ respect to sr filter
        %tic;
        spvals = [];
        if nSpikeWts > 0
           for t = 1:stimLen
                tdiff = t - spindx(spindx < t & spindx >= (t-nSpikeWts));
                for m = 1:nSpikeWts
                  %add indicies and value to sparse matrix values              
                  sval = sum(tdiff == m);
                  if sval > 0
                    spvals = [spvals; m, t, sval];
                  end
                end	    
            end
        end
        %tspike1 = toc;

        %% compute spike-driven part of gradient the same way as stimulus part
        %tic;
        if nSpikeWts > 0
            spMat = sparse(spvals(:, 1), spvals(:, 2), spvals(:, 3), nSpikeWts, stimLen);
            g1Add = (spMat(:, spindx) * rhs');
            g1Spike = g1Spike + g1Add';

            g2Add = (spMat * drate');
            g2Spike = g2Spike + g2Add';        
        end
        %tspike2 = toc;

        %% compute gradient with respect to bias term
        if ~modelParams.ignoreBias           
            g1Bias = g1Bias + sum(irate(spindx) .* drate(spindx));
            g2Bias = g2Bias + sum(drate);            
        end

    end
    
    if ~modelParams.ignoreBias
        g1Stim = [g1Stim g1Bias];
        g2Stim = [g2Stim g2Bias];
    end
    
    g1 = [g1Stim g1Spike];
    g2 = [g2Stim g2Spike]*dt; %g2 is an integral, make sure it's multiplied by dt
    
    %if isfield(modelParams, 'normC')
    %    g1 = g1 / modelParams.normC;
    %    g2 = g2 / modelParams.normC;
    %end
    
    %clear stimulus history matrix early (it can be quite big)
    clear shist;
    
    %% compute derivative for regularizer
    gr = zeros(1, nStimWts + nSpikeWts);
    if modelParams.regularize
        lambda = modelParams.lambda;        
        [modelParams, w] = glppPak(modelParams);
        gr1 = (1 - lambda)*sign(w);
        gr2 = 2*lambda*w;
        gr = modelParams.normC*(gr1 + gr2);
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
    