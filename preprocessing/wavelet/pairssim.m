function [simat,im]=pairssim(cCell,cmpLev)
%function [simat,im]=pairssim(cCell,cmpLev)
%
% Compute the pairwise similarity between all elements of
% [cCell], which are CWT coeff structures. The scale level
% to be used for comparison is specified by [cmpLev].
%
% INPUT:
%  [cCell] = Cell array of CWT coefficient structure that
%            contain 2 fields: [.coef] and [.siz]
% [cmpLev] = vector of wavelet scale levels to be compared
% OUTPUT:
%  [simat] = upper triangular matrix containing similarity 
%            measure between all pairwise elements of [cCell].
%     [im] = cell array of image patches, reconstructed from
%            the corresponding [cCell]
%
% SEE ALSO: cwtssim, dispPairCmp
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Initialization
%--------------------
nC=length(cCell);
cSiz=size(cCell{1},1);
if (nC>20) & (nargout<1)
  warning('pairssim >< Display > 20 patches may be slow');
end
nLev=length(cmpLev);
simat=zeros(nC);



% Compute all pairwise similarity metric
%--------------------
for ii=1:nC
  for jj=ii+1:nC
    simat(ii,jj)=cwtssim(cCell{ii},cCell{jj},cmpLev);
  end  % for jj
end  % for ii



% Compute Inverse Dual Tree CWT for all coefficients
%--------------------
if nargout<1 | nargout>1
  for ii=1:nC
    im{ii}=cwtInv(cCell{ii});
  end  % for ii
end



% Display similarity measures
%--------------------
if nargout<1
  dispPairCmp(im,simat)
end
  
