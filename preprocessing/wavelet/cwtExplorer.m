function [cwtC,im]=cwtExplorer(pix,nLev)
%function [cwtC,im]=cwtExplorer(pix,nLev)
%
% Complex Wavelet Transform (CWT) Explorer. Allow user to click on
% the CWT subimages to turn on=1, or off=0, a wavelet coefficient
% for the blank image [im] or an image the user input [pix].
%      *** Press 'r' or 'R' to reconstruct.
%      *** Press 'q' or 'Q' to quit.
%
% INPUT:
%  [pix] = # of pixel for the sides of the blank image. Or an image
%          in the form of an array of pixel values. (recommend <=32)
% [nLev] = # of wavelet scale to use for decomposition. (recommend <=4)
% OUTPUT:
% [cwtC] = the modified subimages by user's clicking.
%   [im] = image reconstructed by turning on/off certain wavelets.
%
% SEE ALSO: cwtDT, cwtInv, cwtCoef2SubIm, cwtDispSubIm
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input
%--------------------
pSiz=size(pix);
pixIsImg=false;
if all(pSiz>1);  % if pix is an image
  pixIsImg=true;
  pic=pix;
  pix=min(pSiz);
else
  pSiz=[pix,pix];
end

if pix>32
  warning('cwtExplorer >< Larger then 32x32 image may make this explorer run slow');
end
if nLev>4
  warning('cwtExplorer >< [nLve]>4 may make this explorer run slow');
end


% Initialization & setup figures
%--------------------
picZ=zeros(pSiz);
notQuit=true;
hfw=figure;
hfr=figure;
warning off;


% Set up initial subimages
%--------------------
if pixIsImg
  subIc=cwtDT(pic,nLev);   % Working coeff sub image
  subIo=subIc;             % Original coeff sub image
  subIz=cwtDT(picZ,nLev);  % Mask to be flip between 0 & 1
  subIz=cwtZeroLev(subIz,1:nLev+1,1);  % set mask to 1, so when clicked 1st time, it's turn to 0;
else
  subIc=cwtDT(picZ,nLev);
  subIo=cwtZeroLev(subIc,1:nLev+1,1);  % set all original coeff to 1
  subIz=subIc;
end


% Compute sub image sizes
%--------------------
siz=zeros(nLev+1,2);
for ii=1:nLev
  siz(ii,:)=size(subIc{ii}.p75);
end
lowpassSiz=size(subIc{end});
siz(nLev+1,:)=lowpassSiz;


% Display initial figures of sub images
%--------------------
[hfw,hs]=cwtShowStr(subIc,hfw);
axes(hs(end,1));
imagesc(zeros(lowpassSiz));
colormap gray;
axis off; axis image;
ht=text(lowpassSiz(1)/2+0.5,lowpassSiz(2)/2+0.5,['q/r'], ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold', ...
        'FontSize',10);
hs(end,1)=nan;


% Loop for updating CWT Coefficients by clicking on subimages
%--------------------
while notQuit
  warning off;
  % Get Mouse Click / coordinate to change coefficients
  figure(hfw);
  [y,x,b]=ginput(1);
  ha=gca;
  [levIdx,oriIdx]=find(hs==ha);
  x=stickyBound(round(x),1,siz(levIdx,1));
  y=stickyBound(round(y),1,siz(levIdx,2));

  % Update CWT Coefficients on SubImages
  if ~isempty(hs(levIdx,oriIdx)) & b<=3
    axes(hs(levIdx,oriIdx));
    if levIdx<=nLev
      % Get oriented coefficients
      switch oriIdx
        case 1
          subIz{levIdx}.p75(x,y)=mod(subIz{levIdx}.p75(x,y)+1,2);
          subIc{levIdx}.p75(x,y)=subIz{levIdx}.p75(x,y)* ...
                                 subIo{levIdx}.p75(x,y);
          cimage5(subIc{levIdx}.p75);
          axis off; axis image;
        case 2
          subIz{levIdx}.p45(x,y)=mod(subIz{levIdx}.p45(x,y)+1,2);
          subIc{levIdx}.p45(x,y)=subIz{levIdx}.p45(x,y)* ...
                                 subIo{levIdx}.p45(x,y);
          cimage5(subIc{levIdx}.p45);
          axis off; axis image;
        case 3
          subIz{levIdx}.p15(x,y)=mod(subIz{levIdx}.p15(x,y)+1,2);
          subIc{levIdx}.p15(x,y)=subIz{levIdx}.p15(x,y)* ...
                                 subIo{levIdx}.p15(x,y);
          cimage5(subIc{levIdx}.p15);
          axis off; axis image;
        case 4
          subIz{levIdx}.n15(x,y)=mod(subIz{levIdx}.n15(x,y)+1,2);
          subIc{levIdx}.n15(x,y)=subIz{levIdx}.n15(x,y)* ...
                                 subIo{levIdx}.n15(x,y);
          cimage5(subIc{levIdx}.n15);
          axis off; axis image;
        case 5
          subIz{levIdx}.n45(x,y)=mod(subIz{levIdx}.n45(x,y)+1,2);
          subIc{levIdx}.n45(x,y)=subIz{levIdx}.n45(x,y)* ...
                                 subIo{levIdx}.n45(x,y);
          cimage5(subIc{levIdx}.n45);
          axis off; axis image;
        case 6
          subIz{levIdx}.n75(x,y)=mod(subIz{levIdx}.n75(x,y)+1,2);
          subIc{levIdx}.n75(x,y)=subIz{levIdx}.n75(x,y)* ...
                                 subIo{levIdx}.n75(x,y);
          cimage5(subIc{levIdx}.n75);
          axis off; axis image;
      end  % switch oriIdx
      
    else  % if levIdx>nLev
      % Get low pass coefficients
      subIz{levIdx}(x,y)=mod(subIz{levIdx}(x,y)+1,2);
      subIc{levIdx}(x,y)=subIz{levIdx}(x,y)* ...
                         subIo{levIdx}(x,y);
      imagesc(subIc{levIdx});
      colormap gray;
      axis off; axis image;

    end  % if levIdx
  end  % if ~isempty
  
  % Reconstruction from modified CWT coefficients when 'r' or 'R' is press
  if b==114 | b==82
    figure(hfr);  
    im=cwtInv(subIc);
    imagesc(im);
    colormap gray;
    axis off; axis image;
  end

  % Exit Condition = key press 'q' or 'Q'
  if b==113 | b==81
    notQuit=false;
    if nargout>1
      im=cwtInv(cwtC);
    end
    close all;
    warning on;
  end
end  % while


