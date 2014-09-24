function Y = colfilter(X, h)
% function Y = colfilter(X, h)
% Filter the columns of image X using filter vector h, without decimation.
% If length(h) is odd, each output sample is aligned with each input sample
% and Y is the same size as X.
% If length(h) is even, each output sample is aligned with the mid point
% of each pair of input samples, and size(Y) = size(X) + [1 0]; 
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000

[r,c] = size(X);
m = length(h);
m2 = fix(m/2);

% Symmetrically extend with repeat of end samples.
%xe = reflect([(1-m2):(r+m2)], 0.5, r+0.5); % Use 'reflect' so r < m2
% works OK.
% The following is more appropriate for odd-length filters
% JKR 4/26/2001
xe = reflect([(1-m2):(r+m2)], 1, r);

% Perform filtering on the columns of the extended matrix X(xe,:), keeping
% only the 'valid' output samples, so Y is the same size as X if m is odd. 
Y = conv2(X(xe,:),h(:),'valid'); 

return;
