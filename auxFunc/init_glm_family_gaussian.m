function expFam = init_glm_family_gaussian(var)

    if nargin == 0
        var = 1;
    end

    expFam = struct;
    expFam.type = 'gaussian';
    expFam.var = var;
    expFam.meanVar = @(u) var;
    expFam.meanVarInv = @(u) 1 / var;
    expFam.cumFunc = @(x) x.^2 / 2;
    expFam.canFunc = @(u) u;