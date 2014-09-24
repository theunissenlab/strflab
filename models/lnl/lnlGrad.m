function [modelParams, g] = lnlGrad(modelParams, datIdx)
    
    error('Cannot take the gradient of combined linear-nonlinear model...\n');
    
    