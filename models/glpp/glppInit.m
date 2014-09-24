function modelParams = glppInit(numStimChannels, delays, nlType, spikeResponseDelays, sampleRate, ignoreBias)

   if nargin < 6
     ignoreBias = 1;
   end

    %% set basic parameters
    modelParams = struct;
    modelParams.type = 'glpp';
    
    if nargin < 3
        nlType = 'exponential';
    end
    
    if nargin < 4
        spikeResponseDelays = [];
    end
    
    if nargin < 5
        sampleRate = 1000;
    end

    modelParams.ignoreBias = ignoreBias;
    
    modelParams.useLowMemory = 0;
    modelParams.nIn = numStimChannels;
    modelParams.delays = delays;
    modelParams.nWts = length(delays)*numStimChannels + length(spikeResponseDelays) + ~modelParams.ignoreBias;
    
    if strcmp(nlType, 'exponential')
        
        modelParams.nlType = nlType;
        modelParams.nlFunc = @(x) exp(x);   
        modelParams.nlFuncDerivative = @(x) exp(x);

    elseif (strfind(nlType, 'logexp') == 1) & (length(nlType) > 6)
            
        pwr = str2double(nlType(7));

        modelParams.nlType = nlType;
        modelParams.nlFunc = @(x) (log(1 + exp(x))).^pwr;   
        modelParams.nlFuncDerivative = @(x) (pwr*log(1+exp(x)).^(pwr-1)) .* (exp(x) ./ (1 + exp(x)));

    elseif strcmp(nlType, 'harris')
        
        modelParams.nlType = nlType;
        modelParams.nlFunc = @(x) harrisNL(x);
        modelParams.nlFuncDerivative = @(x) dharrisNL(x);
        
    elseif strcmp(nlType, 'custom')        
        modelParams.nlType = 'custom';
        fprintf('[glppInit] make sure you set nlFunc and nlFuncDerivative for custom nonlinearity.\n');
    else
        error('Unknown nonlinearity %s', nlType);
    end
        
    %modelParams.w1 = 1e-6*randn(numStimChannels, length(delays));
    modelParams.w1 = zeros(numStimChannels, length(delays));
    modelParams.b1 = 0;
    modelParams.regularize = 0;
    modelParams.lambda = 0.5;
    modelParams.maxGradMem = 200;
    
    modelParams.l2weight = 1;
    
    modelParams.sampleRate = sampleRate;
    
    spikeResponseWeights = zeros(1, length(spikeResponseDelays));
    modelParams.spikeResponseDelays = spikeResponseDelays;
    modelParams.spikeResponseWeights = spikeResponseWeights;
    
    modelParams.fixedStrf = 0;
    modelParams.useFDGrad = 0;
    modelParams.checkGrad = 0;
    
    modelParams.debugL1 = 0;
    modelParams.debugL2 = 0;
    