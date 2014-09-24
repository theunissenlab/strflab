function [modelParams, w] = glppPak(modelParams)

    w = modelParams.w1(:)';

    if ~modelParams.ignoreBias
        w = [w modelParams.b1];
    end
    
    if ~isempty(modelParams.spikeResponseWeights)
        w = [w modelParams.spikeResponseWeights];
    end
    