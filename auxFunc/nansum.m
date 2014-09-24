function out = nansum(x, varargin)
%NANSUM Sum, ignoring NaNs.
%   M = NANSUM(X) returns the sum of X, treating NaNs as missing
%   values.  For vector input, M is the sum of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the sum of
%   non-NaN elements in each column.  For N-D arrays, NANSUM operates
%   along the first non-singleton dimension.
%
%   NANSUM(X,DIM) takes the sum along dimension DIM of X.	

if nargin == 1
    dim = find(size(x)>1);
else
    dim = varargin{1};
end

nanidx = isnan(x);
x(nanidx)=0;
out = sum(x,dim(1));
