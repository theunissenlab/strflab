function yvals = gompertz(a, b, c, xvals)

    yvals = a*exp(b*exp(c*xvals));
