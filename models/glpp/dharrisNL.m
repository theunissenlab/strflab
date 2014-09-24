function drate = dharrisNL(x)

    zindx = (x <= 0);
    nzindx = ~zindx;
    
    drate = zeros(size(x));
    drate(zindx) = exp(x(zindx));
    drate(nzindx) = ones(size(x(nzindx)));
    
    