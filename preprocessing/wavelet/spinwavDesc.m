function [desc,desParm,desCoef]=spinwavDesc(imSeq,nDesLev)
%function [desc,desParm,desCoef]=spinwavDesc(imSeq,nDesLev)
%
% Compute the complex wavelet descriptor for a sequence
% of image patches. The CW descriptor consist of wavelet
% coefficients to the highest and 2nd highest scale (ie
% coeffs for the largest and 2nd largest wavelets.
%
% INPUT:
%   [imSeq] - sequence of patches. A 3D array of NxNxT.
% OUTPUT:
%    [desc] - Scale-position invariant wavelet descriptor.
% [desParm] - Structure of parameters for the descriptor
%             transform. Needed for descriptor inversion.
%    .imSiz - size of the image feature patch
%     .nLev - # of levels
%  .descIdx - Index of CWT coeff struct used for descriptor
% [desCoef] - cell array of sparse CWT coefficient strct
%
% SEE ALSO: spinwavInv, patchExt, patchNorm
% 
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input & Define constant
% --------------------
if nargin<2
  nDesLev=3;
elseif nDesLev<2
  nDesLev=2;
  warning('spinwavDesc >< # of levels [nDesLev] set to minimum of 2');
end
seqSiz=size(imSeq);
seqLen=seqSiz(3);
nlev=ceil(log2(min(seqSiz(1:2))));


% Make random image for CWT
% --------------------
loopCont=true;
while loopCont
  randIm=randn(seqSiz(1:2));
  rcoef=cwtDT(randIm,nlev);
  loopCont=~isempty(find(rcoef.coef==0));
end


% Determine CWT levels & coeff index to be used as descriptor.
% --------------------
allLev=1:nlev+1;
descLev=allLev(end-nDesLev+1:end);
rmLev=setdiff(allLev,descLev);
zcoef=cwtZeroLev(rcoef,rmLev);
descIdx=find(zcoef.coef~=0);
descLen=length(descIdx);
desc=zeros(descLen,seqLen);


% Compute descriptor
% --------------------
if nargout<=1
  for ii=1:seqLen
    cwtC=cwtDT(imSeq(:,:,ii),nlev);
    desc(:,ii)=cwtC.coef(descIdx);
  end
else  % nargout>1
  desCoef=cell(1,seqLen);
  zeroIdx=find(zcoef.coef==0);
  
  for ii=1:seqLen
    desCoef{ii}=cwtDT(imSeq(:,:,ii),nlev);
    desc(:,ii)=desCoef{ii}.coef(descIdx);
    desCoef{ii}.coef(zeroIdx)=0;
    desCoef{ii}.coef=sparse(desCoef{ii}.coef);
  end
  
  desParm.imSiz=seqSiz(1:2);
  desParm.nLev=nlev;
  desParm.descIdx=descIdx;
end  % if nargout


