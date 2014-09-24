function [modelParams, w] = glmPak(modelParams)

    w = [modelParams.w1(:)' modelParams.b1];
    