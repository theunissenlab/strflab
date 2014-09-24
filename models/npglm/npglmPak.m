function [modelParams, w] = npglmPak(modelParams)

    w = [modelParams.w1(:)' modelParams.b1 modelParams.w2];
    