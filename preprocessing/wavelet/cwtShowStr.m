function [h,hs]=cwtShowStr(cfStr,h)
%function [h,hs]=cwtShowStr(cfStr,h)
%
% Display the extracted wavelet coefficients as sub-image 
% structure array at each scale level.
%
% INPUT:
% [cfStr] = cell array of coefficients at each scale level.
%           Orientation of each level is stored as a structure.
%           So that 
%             [cfStr{N}.p15] = Nth scale level at +15 deg.
%             [cfStr{M}.n45] = Mth scale level at -45 deg.
%           The highest level is the low pass level and has
%           no orientation structure.
%             [cfStr{end}] = low pass filter level. 
% OUTPUT:
%     [h] = figure handle for where to diplay the coefficients.
%    [hs] = axes handle for each subimage of CWT coefficients.
%
% SEE ALSO: cwtDT, cwtInv, cwtZeroLev, cwtZeroOri, cwtZeroBlock
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check input
%--------------------
if nargin<2
  h=figure;
end  % if nargin
nLev=length(cfStr)-1;
hs=nan(nLev+1,6);
siz=zeros(nLev,2);


% Display positive deg CWT coefficients
%--------------------
warning off;
for ii=1:nLev
  figure(h);
  siz(ii,1:2)=size(cfStr{ii}.p15);
  
  % Compute subplots locations for positive orientations
  %--------------------
  subDim=2^ii;
  subDim2=2*subDim;
  totGrid=subDim*subDim;
  subStart=(totGrid-2*subDim)+1;
  subFinal=subStart+subDim+1;
  
  % Display CWT +15, +45 & +75 oriented coefficients
  %--------------------
  hs(ii,1)=tightsubplot(subDim2,subDim,subStart);
  cimage5(cfStr{ii}.p75);
  axis off; axis image;
  
  hs(ii,2)=tightsubplot(subDim2,subDim,subStart+1);
  cimage5(cfStr{ii}.p45);
  axis off; axis image;
  
  hs(ii,3)=tightsubplot(subDim2,subDim,subFinal);
  cimage5(cfStr{ii}.p15);
  axis off; axis image;
  
  % Compute subplots locations for negative orientations
  %--------------------
  subStart=subDim*subDim+2;
  subNext=subStart+subDim;
  
  % Display CWT -15, -45 & -75 deg coefficients
  %--------------------
  hs(ii,4)=tightsubplot(subDim2,subDim,subStart);
  cimage5(cfStr{ii}.n15);
  axis off; axis image;
  
  hs(ii,6)=tightsubplot(subDim2,subDim,subNext-1);
  cimage5(cfStr{ii}.n75);
  axis off; axis image;
  
  hs(ii,5)=tightsubplot(subDim2,subDim,subNext);
  cimage5(cfStr{ii}.n45);
  axis off; axis image;
end


% Display low-pass scale coefficients 
%--------------------
figure(h);
ii=nLev+1;

hs(ii,1)=tightsubplot(subDim2,subDim,subStart-1);
imagesc(cfStr{ii});
colormap gray;
axis off; axis image;

hs(ii,2)=tightsubplot(subDim2,subDim,subStart-subDim-1);
imagesc(cfStr{ii});
colormap gray;
axis off; axis image;


% Set axis display properties
%--------------------
nRow=siz(1,1);
mCol=siz(1,2);
mulFact=round(min(800./[2*nRow,mCol]));
pos=get(h,'Position');
set(h,'Position',[pos(1:2),mulFact*mCol,2*mulFact*nRow]);
warning on;

