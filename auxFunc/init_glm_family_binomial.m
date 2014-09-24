function expFam = init_glm_family_binomial(numTrials)

    expFam = struct;
    expFam.type = 'binomial';
    expFam.numTrials = numTrials;
    expFam.var = 1;
    expFam.meanVar = @(u) u.*(1 - (u ./ numTrials)); 
    expFam.meanVarInv = @(u) (u.*(1 - (u ./ numTrials))).^-1;
    expFam.cumFunc = @(x) numTrials*log(1+exp(x));
    expFam.canFunc = @(u) igfb_canFunc(u, numTrials);

end

function x = igfb_canFunc(u, numTrials)    
    u(u <= 0) = 1e-12;
    x = log(u./(numTrials-u));
end
