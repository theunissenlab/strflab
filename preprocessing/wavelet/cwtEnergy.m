function [eMap,eScale]=cwtEnergy(subIm,gvar,pparm)
%function [eMap,eScale]=cwtEnergy(subIm,gvar,pparm)
%
% Compute the accumulated energy map of the CWT coeff
% of an image.
%
% INPUT:
% [subIm] = CWT wavelet coefficient subimages from cwtCoef2SubIm.
%  [gvar] = variance of the spherical gaussian filter that is
%           used in Gaussian Interpolation of the decimated
%           energy map. Default [gvar]=1.
% [pparm] = Power Parameter that is used to weight the product
%           of the magnitude wavelet coefficients. Default
%           [pparm]=0.25. Reducing this parameter makes the 
%           energy map more sensitive to fine scale structures,
%           but are also more sensitive to noise.
% OUTPUT:
%   [eMap] = Accumulated energy map of CWT coeffients. Obtained
%            by summing [eScale] over the scales.
% [eScale] = Gaussian interpolated energy map at each scale.
%
% SEE ALSO: cwtCoef2SubIm, cwtKeyScale
%
% By Michael Wu  --  waftingpetal@yahoo.com (May 2007)
%
% ====================


% Check Input
%--------------------
if nargin<2
  gvar=1;
end
if nargin<3
  pparm=0.25;
end


% Initialize
%--------------------
nlev=length(subIm)-1;
imDim=2*size(subIm{1}.p15);
imCenter=imDim/2;
fname=fieldnames(subIm{1});
nband=length(fname);

eScale=zeros([imDim,nlev]);
eMap=zeros(imDim);

zeroIm=zeros(imDim);
[gx,gy]=meshgrid([1:imDim(1)]-imCenter(1)-1,[1:imDim(2)]-imCenter(2)-1);


% Compute Scaled Energy Map
%--------------------
for ss=1:nlev
  
  bandDim=size(subIm{ss}.p15);
  bandProd=ones(bandDim);
  for ii=1:nband
    suband=getfield(subIm{ss},fname{ii});
    bandProd=bandProd.*abs(suband);
  end
      
  % Compute Scale Energy Map
  %--------------------
  %enrg=(weightParm.^ss).*((bandProd).^pparm);
  enrg=((bandProd).^pparm);
  
  % Gaussian kernel Fourier interpolation
  %--------------------
  ftE=fftshift(fft2(enrg));
  ftEdim=size(ftE);
  insIdx=imCenter-ftEdim/2;
  ftEpad=zeroIm;
  ftEpad(insIdx(1)+1:insIdx(1)+ftEdim(1),insIdx(2)+1:insIdx(2)+ftEdim(2))=ftE;

  % Gaussian filtering kernel
  %--------------------
  Gz=exp(-0.5*(gx.^2+gy.^2)./(gvar*2.^ss))./(2*pi*gvar*2.^ss);
  ftG=abs(fftshift(fft2(Gz)));
  
  % Inverse Fourier Transform
  %--------------------
  eScale(:,:,ss)=abs(ifft2(ifftshift(ftEpad.*ftG)));
end  % for ii


% Compute Accumulated Energy Map
%--------------------
eMap=sum(eScale,3);

