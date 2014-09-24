function [oriCoef]=crvGetOri(ccf,oriDeg)
%function [oriCoef]=crvGetOri(ccf,oriDeg)
% 
% Extract the curvelet coefficients at an orientation [oriDeg] in
% degrees. If no output is assinged, it will plot the real part of
% the coefficients.
%
% INPUT:
% [ccf]     : cell array of curvelet coefficients at each scale [s],
%             orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
% [oriDeg]  : orientation we wish to extract at different scales.
%             Horizontal = 0 deg, and vertical = 90 deg. 
%
% OUTPUT:
% [oriCoef] : cell array of curvelet coefficients at all scale with
%             orientation [oriDeg]. oriCoef{w} = coefficient at scale 
%             [scaleLev] and orientation [w].
%
% SEE ALSO: crvShowOri, crvGetLev, crvGetBlock, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
oriDeg=mod(oriDeg,360);
nLev=length(ccf);
finestLen=length(ccf{nLev});
if finestLen>1
  oriCoef=cell(1,nLev);  % finest level is curveleet
else
  oriCoef=cell(1,nLev-1);  % finest level is wavelet
end
nLevwOri=length(oriCoef);


% Extract coefficients at orientation [oriDeg]
%--------------------
oriRng=nan(nLevwOri,2);
for ii=2:nLevwOri
  nOri=length(ccf{ii});
  % Convert Orientation in deg to curvelet orientation index
  % oriIdx may contain 2 index if desired orientation lies at the
  % boundary of two index. But only the first will be use.
  [oriIdx,oriRng(ii,:)]=crvOridx(nOri,oriDeg);
  oriCoef{ii}=ccf{ii}{oriIdx{1}};
end


% Plot coefficients when no output is assigned.
%--------------------
if nargout<1
  hFig=figure;
  crvShowOri(oriCoef,oriRng,hFig);
end

