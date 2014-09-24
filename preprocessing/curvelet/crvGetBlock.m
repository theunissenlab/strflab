function [blkCoef]=crvGetBlock(ccf,levScale,oriDeg) 
%function [blkCoef]=crvGetBlock(ccf,levScale,oriDeg) 
% 
% Extract the curvelet coefficients at scale [levScale] and orientation
% [oriDeg] in degrees. If no output is assinged, it plots the real part of
% the coefficient.
%
% INPUT:
%      [ccf] : cell array of curvelet coefficients at each scale [s],
%              orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
% [levScale] : The scale to extract the coefficient. 
%              1 = coarest scale, increasing [scaleLev] = finer scale.
%   [oriDeg] : orientation we wish to extract at different scales.
%              Horizontal = 0 deg, and vertical = 90 deg. 
% OUTPUT:
%  [blkCoef] : matrix of curvelet coefficients extracted.
%
% SEE ALSO: crvGetLev, crvGetOri, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Check Scale Info
%====================
nLev=length(ccf);
finestLen=length(ccf{nLev});


% Extract coef @ scale [scaleLev]
%====================
if levScale<=nLev
  levCoef=ccf{levScale};
else
  error('crvGetBlock >< [levScale] must be <= # of levels in the [ccf]');
end


% Check Orientation Info
%====================
levCoefLen=length(levCoef);
oriDeg=mod(oriDeg,360);
oriIdx=crvOridx(levCoefLen,oriDeg);


% Extract coef with orientation [oriDeg]
%====================
blkCoef=levCoef{oriIdx{1}};
if levCoefLen<=1 & nargin>2
  warning('crvGetBlock >< no oriented coefficients at DC & finest scale');
end


% Plot coefficients when no output is assigned.
%====================
if nargout<1
  hFig=gcf; 
  imagesc(real(blkCoef));
  axis image; colormap gray; axis off;
end

