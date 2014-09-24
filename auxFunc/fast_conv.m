function a = fast_conv(stim,filter,d,T,S)
%function a = fast_conv(stim,filter,d,T,S)
%
%  Does a fast convolution between the T by S stimulus and the S by d filter
%
% INPUT:
%   [stim] = stimulus vector or matrix
% [filter] = filter to be convolved with the stimulus
%      [d] = # of columns of filter
%      [T] = # of rows of stimulus
%      [S] = # of columns of stimulus & # of rows of filter
%
% OUTPUT:
%      [a] = stimulus convolved with the filter

        prod = stim*filter;
        out = zeros(d,T);
        for jj = 1:d
            out(:,jj) = [prod((1+(jj-1)*T):(1-T):1) zeros(1,d-jj)];
        end

        for jj = 0:(T-d)
            out(:,jj + d) = prod((T*(d-1) + jj + 1):(1-T):(d+jj));
        end
        a = sum(out);
