function [modelParams, w] = leglmPak(modelParams)

    w = [modelParams.w1(:)' modelParams.b1 modelParams.m];
    