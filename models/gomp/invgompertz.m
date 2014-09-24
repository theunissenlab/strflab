function xvals = invgompertz(a, b, c, yvals)

    aeps = 1e-8;
    zeps = 1e-8;
    yvals(yvals >= a) = a - aeps;
    yvals(yvals == 0) = zeps;
    xvals = log( log(yvals / a) / b) / c;
