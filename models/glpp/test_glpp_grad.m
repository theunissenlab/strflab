function test_glpp_grad()

    

    sampleRate = 1000;
    duration = 20;
    nDataPoints = round(duration * sampleRate);

    nStimDelays = 3;
    numStimChannels = 2;
    stimDelays = 0:(nStimDelays-1);
    nSpikeDelays = 5;
    spikeDelays = 0:(nSpikeDelays-1);
    
    nlType = 'exponential';
        
    %% initialize model
    ignoreBias = 0;
    modelParams = glppInit(numStimChannels, stimDelays, nlType, ...
                		       spikeDelays, sampleRate, ignoreBias);

    w = randn(numStimChannels, nStimDelays);
    w = (w / norm(w(:))) * 4;
    b = log(0.02);
    
    modelParams.w1 = w;
    modelParams.b1 = b;
    modelParams.normC = 8000;
    modelParams.spikeResponseWeights = randn(1, nSpikeDelays);
    modelParams.regularize = 1;
    modelParams.lambda = 0.5;
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels);
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);    

    numRealTrials = 20;    
    datIdx = 1:nDataPoints;
    [modelParams, modelResponseReal, fullResponse] = glppFwd(modelParams, datIdx, numRealTrials);

    sint = 1/sampleRate;
    t = 0:sint:(nDataPoints-1)*sint;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    
    spikeTimesCell = cell(numRealTrials, 1);
    for k = 1:numRealTrials
      strials = fullResponse.spikeTrials(k, :);
      stimes = (find(strials > 0) - 1) / sampleRate;
      spikeTimesCell{k} = stimes;      
    end
    strfData(stim, spikeTimesCell, groupIndex);
    
    datIdx = find(groupIndex == 1);
    
    [modelParams, g] = glppGrad(modelParams, datIdx);    
    [modelParams, gfd] = glppGradFD(modelParams, datIdx);
    
    g
    gfd
    gdiff = norm(g - gfd)
    
    
    

    