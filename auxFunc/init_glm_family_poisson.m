function expFam = init_glm_family_poisson()

    expFam = struct;
    expFam.type = 'poisson';
    expFam.var = 1;
    expFam.meanVar = @(u) u;
    expFam.meanVarInv = @(u) 1 ./ u;
    expFam.cumFunc = @(x) exp(x);
    expFam.canFunc = @igfp_canFunc;
    
    
end

function x = igfp_canFunc(u)
    u(u <= 0) = 1e-12;
    x = log(u);
end
