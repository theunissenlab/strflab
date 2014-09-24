function [dss]=crvPlotSS(ssRep,hFig)
%function hsp=crvPlotSS(ssRep,hFig)
% 
% Plot the image at different scale space representation.
%
% INPUT:
%   ssRep : cell array of images at different curvelet scale space
%           representation. Each cell is a image filtered at different
%           scale. ssRep{1} most coarse scale.
%
% OUTPUT:
%   hsp   : subplot handles for the plots.
%
% SEE ALSO: crvSSrep, crvScale, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%====================
if nargin>1
    figure(hFig);
else
    hFig=figure;
end
nScale=length(ssRep);


% Determine subplot Grid Dimension
%====================
spCol=ceil(sqrt(nScale));
if nScale>(spCol*(spCol-1))
    spRow=spCol;
else
    spRow=spCol-1;
end


% Plot Scale Space Representation
%====================
hsp=zeros(nScale,1);
for ss=1:nScale
    figure(hFig);
    hsp(ss)=subplot(spRow,spCol,ss);
    imagesc(ssRep{ss});
    axis off; axis image; colormap gray;    
end


% Plot scale space derivatives (difference)
%====================
hFig2=figure;
hsd=zeros(nScale-1,1);
dss=cell(1,nScale-1);
for ss=1:nScale-1
    figure(hFig2);
    dss{ss}=ssRep{ss}-ssRep{ss+1};
    hsp(ss)=subplot(spRow,spCol,ss);
    imagesc(dss{ss});
    axis off; axis image; colormap gray;    
end


