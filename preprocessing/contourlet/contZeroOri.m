function [cf]=contZeroOri(cf,oriDeg,c)
%function [cf]=contZeroOri(cf,oriDeg,c)
% 
% Remove the contourlet coefficients at orientations [oriDeg] in 
% degrees.
%
% INPUT:
%     [cf] : cell array of curvelet coefficients at each scale [s],
%            orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
% [oriDeg] : orientations we wish to remove (set to 0) at different 
%            scales. Horizontal = 0 deg, and vertical = 90 deg. 
%      [c] : a scalar, which the coeff in [levScale] are set to.
%            Default=0.
% OUTPUT:
%     [cf] : same as above, but all coeff at orientations [oriDeg]
%            are set to zero.
%
% SEE ALSO: 
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
oriDeg=mod(oriDeg,360);
rmOri=unique(oriDeg);
nRmOri=length(rmOri);

nLev=length(ccf);
finestLen=length(ccf{nLev});
if finestLen==1
    nLev=nLev-1;
end


% Remove coefficients at orientation [oriDeg]
%--------------------
for ss=2:nLev
  nOri=length(ccf{ss});
  oriIdx=crvOridx(nOri,oriDeg);
  uniqIdx=unique(cell2mat(oriIdx))';
  for ii=uniqIdx
    ccf{ss}{ii}(:)=0;
  end
end

