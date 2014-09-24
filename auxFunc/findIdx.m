function [idx]=findIdx(idxLis,idxSet)
%function [idx]=findIdx(idxLis,idxSet)
%
% Find the index of multiple elements in [idxLis] of the vector 
% [indSet].
%
% INPUT:
%  idxLis = vector of elements that are members of [idxSet]
%  idxSet = vector of elements where the index is sought.
%
% OUTPUT:
%   [idx] = vector of sorted index of all elements of [idxLis]
%           in [idxSet].
%
% EXAMPLE:
%   g=[1 4 2 1 2 2 3 3 1 4 4 2 4 1 2];
%   findIdx(4,g)
%   >>  2  10  11  13                      % same as find(4==g);
%   findIdx([1 2],g)
%   >>  1   3   4   5   6   9  12  14  15  % kind of like 
%                                          % find([1 2]==g);
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jul 2007)
%
% ====================


idxLen=length(idxLis);
setLen=length(idxSet);

idx=[];
for ii=1:idxLen
  idx=[idx,find(idxLis(ii)==idxSet)];
end
idx=sort(idx);

