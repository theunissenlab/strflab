function [ssRep]=crvSSrep(ccoef,mRow,nCol)
%function [ssRep]=crvSSrep(ccoef,mRow,nCol)
% 
% Compute the curvelet scale space representation of an image by
% reconstructing the image at each scale from the curvelet coefficients at
% the corresponding scales.
%
% INPUT:
%   ccoef : cell array of curvelet coefficients at each scale [s],
%           orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
%   mRow  : height of the image
%   nCol  : width of the image
%
% OUTPUT:
%   ssRep : cell array of images at different curvelet scale space
%           representation. Each cell is a image filtered at different
%           scale. ssRep{1} most coarse scale.
%
% SEE ALSO: crvPlotSS, crvScale, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%====================
nScale=length(ccoef);
if nargin<2
    coefLen=length(ccoef{nScale});
    if coefLen~=1
        error('crvSSrep >< must provide [mRow] & [nCol]');
    else
        coefSize=size(ccoef{end}{end});
        mRow=coefSize(1);
        nCol=coefSize(2);
    end
end


% Init a zero curvelet coefficients
%====================
zeroCoef=cell(1,nScale);
for ss=1:nScale
    oriSize=size(ccoef{ss});
    zeroCoef{ss}=cell(oriSize);
    for oo=1:oriSize(2)
        zeroCoef{ss}{oo}=zeros(size(ccoef{ss}{oo}));
    end
end


% Compute the scale space rep by coefficients at each curvelet scale.
%====================
for ss=1:nScale
    tmpCoef=zeroCoef;
    tmpCoef{ss}=crvScale(ccoef,ss);
    ssRep{ss}=real(crvInv(tmpCoef,mRow,nCol));
end


% Plot scale space rep when no output is assigned.
%====================
if nargout<1
    hFig=figure;
    crvPlotSS(ssRep,hFig);
end


