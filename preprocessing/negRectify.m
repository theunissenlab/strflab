function x = negRectify(x)

x(x>0) = 0;
x = -x;