function [modelParams, invResp] = gompInvert(modelParams, resp)
    
    if ~strcmp('gomp', modelParams.type)
        error('Cannot invert a non-gomp model!');
    end

    a = modelParams.a;
    b = modelParams.b;
    c = modelParams.c;
    
    invResp = invgompertz(a, b, c, resp);
    