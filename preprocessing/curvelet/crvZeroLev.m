function [ccf]=crvZeroLev(ccf,levScale)
%function [ccf]=crvZeroLev(ccf,levScale)
% 
% Remove the curvelet coefficients at scales specified by [levScale].
%
% INPUT:
%      [ccf] : cell array of curvelet coefficients at each scale [s],
%              orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
% [levScale] : The scales at which we zero-out the coefficient. 
%              1 = coarest scale, increasing [levScale] = finer scale.
%
% OUTPUT:
%      [ccf] : same as above, but with the coeff at scale [levScale]
%              set to zero.
%
% SEE ALSO: crvGetOri, crvGetLev, crvGetBlock, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
nLev=length(ccf);
rmLev=unique(levScale);
nRmLev=length(rmLev);
if rmLev(end)>nLev | rmLev(1)<1
  error('crvZeroLev >< values in [scaleLev] must be >=1 & <= # of scales in [ccoef]');
end


% Remove coefficients at scale [scaleLev]
%--------------------
for ss=rmLev
  nOri=length(ccf{ss});
  for ii=1:nOri
    ccf{ss}{ii}(:)=0;
  end  % for ii
end  % for ss


