function modelParams = glmInit(nIn, delays, family, outputNL, dispersion, qval, qwt)

    if nargin < 7
        qwt = 0;
    end
    if nargin < 6
        qval = 1;
    end        
    
    modelParams.type = 'glm';
    
    modelParams.family = family;
    modelParams.nIn = nIn;
    modelParams.delays = delays;
    modelParams.outputNL = outputNL;
    modelParams.dispersion = dispersion;
    modelParams.w1 = zeros(nIn, length(delays));
    modelParams.b1 = 0;
    modelParams.nWts = nIn*length(delays) + 1;
    
    modelParams.qval = qval;
    modelParams.qwt = qwt;
    
    %% check qval
    if qval >= 2 || qval < 1
        error('For regularization to work, we need 1 < qval <= 2');        
    end
    