function rate = harrisNL(x)

    zindx = (x <= 0);
    nzindx = ~zindx;
    
    rate = zeros(size(x));
    rate(zindx) = exp(x(zindx));
    rate(nzindx) = 1 + x(nzindx);