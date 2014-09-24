function [A, idxr, idxc] = getSubMatrix(hsize, locality, sdimension);
%function A = getSubMatrix(hsize, sdimension)
%
% A function to extract the relevant parts of a covariance matrix
%
% INPUT:
%      [hsize] = vector of sizes for each dimension of the stimulus
%   [locality] = determines how far away to consider pixel interaction terms
% [sdimension] = determines which dimensions to get covariace for
%
% OUTPUT:
%          [A] = sparse matrix of valid terms
%       [idxr] = index of the rows of the valid terms of A
%       [idxc] = index of the columns of the valid terms of A
%
% EXAMPLE:
%   A = getSubMatrix([20 20 5], 2, [1 1 0]);
%

hsize = hsize(logical(sdimension));

if length(hsize) <= 4
  hsize = [hsize ones(1,4-length(hsize))];
end

if sdimension(3)
    localityt = locality;
    s = ones(locality, locality, localityt);
    scent = ceil(locality/2);
    ssize = size(s);
else
    localityt = 3;
    s = zeros(locality, locality, localityt);
    s(:,:,2) = ones(locality);
    scentr = ceil(locality/2);
    scent = 2;
    ssize = size(s);
end


if exist('sdimension') 
  for sdir=1:3
      if sdimension(sdir)
          ssubs{sdir} = 1:locality;  
      else
          ssubs{sdir} = scent;  
      end
  end
  s2 = zeros(locality, locality, localityt);
  s2(ssubs{1}, ssubs{2}, ssubs{3}) = s(ssubs{1}, ssubs{2}, ssubs{3});
  s = s2;
end



nonzeroind = find(s);
nonzeroval = s(nonzeroind);
numnonzeros = length(nonzeroind);

for jj=1:numnonzeros
  [sind(jj,1) sind(jj,2) sind(jj,3)] = ind2sub(ssize, nonzeroind(jj));
end

if sdimension(3)
    sind = sind - scent;
else
    sind = sind(:,1:2) - scentr;
    sind(:,3) = 0;
end

dimension3 = prod(hsize(1:3));
dimension4 = prod(hsize);

hsize3 = hsize(1:3);

A = sparse([],[],[],dimension4, dimension4, dimension4*length(find(s))*hsize(4));
d3=1:dimension3;
[d(:,1) d(:,2) d(:,3)] = ind2sub(hsize(1:3), d3);
  for jj=1:numnonzeros
    thisind = bsxfun(@plus,d,sind(jj,:));

    inidx = find(checkwithin(hsize3, thisind)>0);
    ind = sub2ind(hsize, thisind(inidx,1), thisind(inidx,2), thisind(inidx,3));
    indA = sub2ind(size(A), inidx, ind);
    A(indA) = 1;

  end


A = triu(A);

[idxr, idxc] = find(A>0);
idxr = uint32(idxr);
idxc = uint32(idxc);

return;


%%%-----------------------------
function out = checkwithin(hsize, indarray)
    
    out = ones([size(indarray,1) 1]);
    [r,c] = find(indarray<=0);
    out(r) = 0;

    indarray2 = bsxfun(@minus,indarray,hsize);
    [r,c] = find(indarray2>0);
    out(r) = 0;
