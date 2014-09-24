function out = nanmean(x, varargin)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.	

if nargin == 1
    dim = min(find(size(x)>1));
else
    dim = varargin{1};
end

nanidx = isnan(x);
x(nanidx)=0;
dimsize = sum(~nanidx,dim(1));
out = sum(x,dim(1)) ./ dimsize;
