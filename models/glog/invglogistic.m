function x = invglogistic(y, B, M)

    deps = 1e-6;
    y(y == 1) = 1 - deps;
    y(y == 0) = deps;

    x = (-log(y.^-1 - 1) / B) + M;
    