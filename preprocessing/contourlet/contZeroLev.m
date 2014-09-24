function [cf]=contZeroLev(cf,levScale,c)
%function [cf]=contZeroLev(cf,levScale,c)
%
% Remove the contourlet coefficients at scales specified by [levScale].
%
% INPUT:
%       [cf] : cell array of contourlet coefficients at each scale [s],
%              orientation [w], and position [x,y]. cf{s}{w}[x,y].
% [levScale] : The scales at which we zero-out the coefficient. 
%              1 = coarest scale, increasing [levScale] = finer scale.
%        [c] : a scalar, which the coeff in [levScale] are set to.
%              Default=0.
% OUTPUT:
%       [cf] : same as above, but with the coeff at all levels in
%              [levScale] set to zero.
%
% SEE ALSO: 
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Check input
%--------------------
if nargin<3
  c=0;
end


% Init Parameters
%--------------------
nLev=length(cf);
rmLev=unique(levScale);
nRmLev=length(rmLev);
if rmLev(end)>nLev | rmLev(1)<1
  error('contZeroLev >< values in [scaleLev] must be >=1 & <= # of scales in [cf]');
end


% Remove coefficients at scale [scaleLev]
%--------------------
for ss=rmLev
  if iscell(cf{ss})
    nOri=length(cf{ss});
    for ii=1:nOri
      cf{ss}{ii}(:)=c;
    end  % for ii
  else
    cf{ss}(:)=c;
  end
end  % for ss

