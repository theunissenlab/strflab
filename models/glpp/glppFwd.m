function [modelParams, modelResponse, fullResponse] = glppFwd(modelParams, datIdx, numTrials)

    global globDat;
    
    if nargin < 3
        numTrials = 1;
    end
    
    %% compute stimulus convolution
    stim = globDat.stim(datIdx, :);
    
    groupIndex = globDat.groupIdx;
    if ~isempty(datIdx)
        stim = globDat.stim(datIdx, :);

        groupIndex = globDat.groupIdx(datIdx);
        uindx = unique(groupIndex);
        %make sure groups go from 1 to length(uindex) instead of 1,2,5,7,..
        %in the event of held out datasets not included by datIdx
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
            
    if nargout > 2
        fullResponse = struct;
        fullResponse.numTrials = numTrials;
        fullResponse.stimCurrent = stimCurrent;
        fullResponse.spikeCurrent = zeros(numTrials, length(stimCurrent));
        fullResponse.nonlinearResponse = zeros(numTrials, length(stimCurrent));
    end
    
    %% simulate spike trials
    spikeTrials = zeros(numTrials, length(globDat.stim(datIdx)));
    swts = modelParams.spikeResponseWeights;
    f = @(x) modelParams.nlFunc(x);
    dt = 1 / modelParams.sampleRate;
    
    nDatasets = length(unique(groupIndex));
    
    for k = 1:numTrials
        
        sumLinearResponse = zeros(1, length(stimCurrent));
    
        %next spike time is from a homogeneous poisson process w/ rate=1
        %because of time-rescaling
        nextSpikeTime = exprnd(1);
        tOffset = 0;
        
        %generate a spike train for each dataset (concatenated along the way)
        for d = 1:nDatasets
        
            rateSum = 0;
            rng = find(groupIndex == d);
            linResp = stimCurrent(rng);
            respLen = length(linResp);
            
            % stimulate the model at each time point
            for t = 1:length(linResp)
                %keep track of the integral of the rate since last spike,
                %a homogeneous poisson process w/ rate=1
                rateSum = rateSum + f(linResp(t))*dt;
                if rateSum >= nextSpikeTime

                    %reset the sum and generate the next spike time
                    spikeTrials(k, t+tOffset) = 1;
                    rateSum = 0;
                    nextSpikeTime = exprnd(1);

                    %add the spike response current to the linear response
                    %for time points after this
                    if ~isempty(swts)
                        lindx = min(t+length(swts), respLen);
                        scrng = (t+1):lindx;
                        swtsToUse = swts(1:length(scrng));
                        linResp(scrng) = linResp(scrng) + swtsToUse;
                    end
                    
                    %record spike current if requested
                    if nargout > 2 && ~isempty(swts)
                        trng = tOffset + scrng;
                        fullResponse.spikeCurrent(k, trng) = fullResponse.spikeCurrent(k, trng) + swtsToUse;
                    end
                end
            end 

            sumLinearResponse(rng) = linResp;
            
            %increment the dataset offset
            tOffset = tOffset + respLen;
            
            %figure; hold on;
            %plot(fullResponse.spikeCurrent, 'g-');
            
        end
        
        %record nonlinear response
        if nargout > 2            
            fullResponse.nonlinearResponse(k, :) = f(sumLinearResponse);
            %{
            figure; hold on;
            subplot(3, 1, 1);
            plot(sumLinearResponse); axis tight;
            title('sumLinearResponse');
            subplot(3, 1, 2);
            plot(fullResponse.nonlinearResponse(k, :)); axis tight;
            subplot(3, 1, 3);
            hist(fullResponse.nonlinearResponse(k, :)); axis tight;
            title('nonlinear response');
            %}
        end
        
    end

    modelResponse = sum(spikeTrials, 1) / numTrials;
    
    if nargout > 2
        fullResponse.spikeTrials = spikeTrials;
        fullResponse.modelResponse = modelResponse;
    end
