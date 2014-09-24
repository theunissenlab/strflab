function modelParams = gompInit(saturationValue)

    modelParams = struct;
    modelParams.type = 'gomp';
    modelParams.nWts = 2;
    modelParams.a = saturationValue;
    modelParams.b = -1;
    modelParams.c = -1;
    