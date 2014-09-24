function [modelParams, w] = lnlPak(modelParams)

    [linm, linw] = strfPak(modelParams.linModel);
    [nlm, nlw] = strfPak(modelParams.nlModel);
    w = [linw nlw];
    