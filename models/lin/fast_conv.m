function out = fast_conv(stim,filter,d,T)
prod = stim*filter;
out = zeros(d,T);
for jj = 1:d
    out(:,jj) = [prod((1+(jj-1)*T):(1-T):1) zeros(1,d-jj)];
end

for jj = 0:(T-d)
    out(:,jj + d) = prod((T*(d-1) + jj + 1):(1-T):(d+jj));
end
out = sum(out);