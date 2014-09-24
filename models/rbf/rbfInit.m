function modelParams = rbfInit(func, dim, nfuncs)

    modelParams.type = 'rbf';
    modelParams.n = dim;
    modelParams.nfuncs = nfuncs;
    modelParams.centers = zeros(nfuncs, dim);
    modelParams.func = func;
    modelParams.w1 = zeros(1, nfuncs);
    modelParams.nWts = nfuncs;
    