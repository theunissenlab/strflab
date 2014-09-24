function test_glpp1()

    sampleRate = 1000;
    duration = 20;
    nDataPoints = round(duration * sampleRate);

    nStimDelays = 75;
    numStimChannels = 60;
    stimDelays = 0:(nStimDelays-1);
    nSpikeDelays = 50;
    spikeDelays = 0:(nSpikeDelays-1);
    
    nlType = 'exponential';
        
    %% initialize model
    ignoreBias = 0;
    modelParamsReal = glppInit(numStimChannels, stimDelays, nlType, ...
                		       spikeDelays, sampleRate, ignoreBias);

    w = randn(numStimChannels, nStimDelays);
    w = (w / norm(w(:))) * 4;
    b = log(0.02);
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = b;
    modelParamsReal.spikeResponseWeights = randn(1, length(spikeDelays));
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels);
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);    

    numRealTrials = 20;    
    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = glppFwd(modelParamsReal, datIdx, numRealTrials);

    sint = 1/sampleRate;
    t = 0:sint:(nDataPoints-1)*sint;
    
    figure; hold on;
    subplot(4, 1, 1);
    plot(t, stim, 'k-');
    title('Stimulus');
    subplot(4, 1, 2);
    plot(t, fullResponse.stimCurrent, 'b-');
    title('Linear Response');
    subplot(4, 1, 3); hold on;
    exampleSets = [1:10:numRealTrials];
    for k = 1:length(exampleSets)
      plot(t, fullResponse.nonlinearResponse(exampleSets(k), :));
    end
    title('Example Poisson Rates');
    subplot(4, 1, 4); hold on;
    plot(t, modelResponseReal, 'k-');
    title('PSTH');
        
    %% Set options for optimization method
    optOptions = trnSCG;
    %optOptions.minErrChange = 1e-2;
    %optOptions.minStepSize = 1e-2;
    optOptions.display = 1;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    
    spikeTimesCell = cell(numRealTrials, 1);
    for k = 1:numRealTrials
      strials = fullResponse.spikeTrials(k, :);
      stimes = (find(strials > 0) - 1) / sampleRate;
      spikeTimesCell{k} = stimes;      
    end
    strfData(stim, spikeTimesCell, groupIndex);
    
    %% split responses for coherence
    preprocRespParams = struct;         
    preprocRespParams.split = 1;
    preprocRespParams.units = 's';     
    preprocRespParams.binSize = 0.001;  
    preprocRespParams.stimLengths = [duration];
    [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses({spikeTimesCell}, preprocRespParams);
    
    %% Initialize LNP spike-response model
    modelParams = glppInit(numStimChannels, stimDelays, nlType, ...
                           spikeDelays, sampleRate, modelParamsReal.ignoreBias);
    
    fprintf('Using perturbed initial guess...\n');
    modelParams.w1 = randn(size(modelParamsReal.w1))*1e-1 + modelParamsReal.w1;
    mr = modelResponseReal;
    mr(modelResponseReal <= 0) = 1e-3;
    modelParams.b1 = mean(log(mr));
    
    %% figure out a good normalization constant
    [modelParams, err] = glppErr(modelParams, find(groupIndex == 1));
    err
    modelParams.normC = abs(err);
    modelParams.checkGrad = 0;    
    modelParams.useFDGrad = 0;
    
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    %% look at response
    numTrials = 100;
    [modelParamsTrained, predResp] = glppFwd(modelParamsTrained, datIdx, numTrials);
    predResp = (predResp/max(predResp)) * max(modelResponseReal);

    realStrf = modelParamsReal.w1
    predStrf = modelParamsTrained.w1
    strfDiff = norm(realStrf - predStrf);
    figure; hold on;
    imagesc(realStrf - predStrf);
    axis tight; colorbar;
    title(sprintf('STRF Diff | norm=%0.6f', strfDiff));
  
    realw0 = modelParamsReal.b1
    predw0 = modelParamsTrained.b1
    
    realSwts = modelParamsReal.spikeResponseWeights
    predSwts = modelParamsTrained.spikeResponseWeights
    
    figure; hold on;
    plot(modelResponseReal, 'k-', 'LineWidth', 2);
    plot(predResp, 'r-', 'LineWidth', 2);
    axis tight;
    legend('Real', 'Predicted');
    %axis([1 length(modelResponseReal) 0 1]);
    
    
    %% compute coherence    
    infoFreqCutoff = 90;
    infoWindowSize = 0.250;
    [cBound, cModel] = compute_coherence_full(predResp, modelResponseReal, wholeSplitResponse(1, :), wholeSplitResponse(2, :), sampleRate, numRealTrials, infoFreqCutoff, infoWindowSize);
        
    performanceRatio = cModel.info / cBound.info

    %% plot information
    figure; hold on;
    plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
    plot(cBound.f, cBound.cUpper, 'b-', 'LineWidth', 2);
    plot(cBound.f, cBound.cLower, 'r-', 'LineWidth', 2);

    plot(cModel.f, cModel.c, 'k--', 'LineWidth', 2);
    plot(cModel.f, cModel.cUpper, 'b--', 'LineWidth', 2);
    plot(cModel.f, cModel.cLower, 'r--', 'LineWidth', 2);
    theTitle = sprintf('Info=%0.0f out of %0.0f bits | Ratio=%0.2f', ...
                       cModel.info, cBound.info, performanceRatio);
    title(theTitle);
    axis([min(cBound.f), max(cBound.f), 0, 1]);
    