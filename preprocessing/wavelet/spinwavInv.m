function [desMov,desCoef]=spinwavInv(desc,desParm)
%function [desMov,desCoef]=spinwavInv(desc,desParm)
%
% Invert the 2 upper scale wavelet coefficients to visualize
% the spinwavDesc features in image space.
%
% INPUT:
%    [desc] - set of spinwave descriptor DxN. Where N = # of
%             features. D=feature dimensionality.
% [desParm] - spinwav descriptor parameters containing fields
%    .imSiz - size of the image feature patch
%     .nLev - # of levels
%  .descLen - length (dimensionality) of descriptor
%  .descIdx - Index of CWT coeff struct used for descriptor
%  .zeroIdx - Index of CWT coeff struct zero out
% OUTPUT:
%  [desMov] - a 3D array of SxSxN, where S is the size of the
%             image patch. N = # of features patches.
% [desCoef] - cell array of CWT coefficient structures that
%             contains [.coef] and [.siz] fields.
%
% SEE ALSO: spinwavDesc, patchExt, patchNorm
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input
% --------------------
if iscell(desc)
  desc=cell2mat(desc);
end  % if iscell
[descDim,nDesc]=size(desc);


% Initialization
% --------------------
rcoef=cwtDT(randn(desParm.imSiz),desParm.nLev);
rcoef.coef(:)=0;
desMov=zeros(desParm.imSiz(1),desParm.imSiz(2),nDesc);


% Invert the wavelet descriptor
% --------------------
if nargout>1
  desCoef=cell(nDesc,1);
  for ii=1:nDesc
    rcoef.coef(desParm.descIdx)=desc(:,ii);
    desMov(:,:,ii)=cwtInv(rcoef);
    desCoef{ii}=rcoef;
  end  % for ii
else
  for ii=1:nDesc
    rcoef.coef(desParm.descIdx)=desc(:,ii);
    desMov(:,:,ii)=cwtInv(rcoef);
  end  % for ii
end  % if nargout


