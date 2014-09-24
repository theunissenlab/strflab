function [levCoef]=crvGetLev(ccf,levScale)
%function [levCoef]=crvGetLev(ccf,levScale)
% 
% Extract the curvelet coefficients at scale [levScale]. If no 
% output is assinged, it will plot the real part of the coefficients.
%
% INPUT:
%   [ccoef] : cell array of curvelet coefficients at each scale [s],
%             orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
%     [lev] : The scale to extract the coefficient. 
%             1 = coarest scale, increasing [scaleLev] = finer scale.
%
% OUTPUT:
% [levCoef] : cell array of curvelet coefficients at scale [lev].
%             levCoef{w} = coefficient at scale [lev] and orientation 
%             with orientation index [w].
%
% SEE ALSO: crvShowLev, crvGetOri, crvGetBlock, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
nlev=length(ccf);
finestLen=length(ccf{nlev});


% Extract coefficients at scale [scaleLev]
%--------------------
if levScale<=nlev
  levCoef=ccf{levScale};
else
  error('crvGetLev >< [levScale] must be <= total # of levels in [ccf]');
end


% Plot coefficients when no output is assigned.
%--------------------
if nargout<1
  hFig=figure;
  crvShowLev(levCoef,hFig);
end

