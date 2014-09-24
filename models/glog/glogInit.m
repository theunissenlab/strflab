function modelParams = glogInit()

    modelParams = struct;
    modelParams.type = 'glog';
    modelParams.nWts = 2;
    modelParams.M = 0;
    modelParams.B = abs(randn()*1.5);
    