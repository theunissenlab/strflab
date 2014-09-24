function y = glogistic(x, B, M)

    Q = 1;
    v = 1;
    y = (1 + Q*exp(-B*(x-M))) .^ (-1/v);
    
    