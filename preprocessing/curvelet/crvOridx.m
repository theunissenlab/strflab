function [oriIdx,oriRng]=crvOridx(nOri,oriDeg,opt)
%ffunction [oriIdx,oriRng]=crvOridx(nOri,oriDeg,opt)
% 
% Computes the curvelet orientation index that correspond to a given
% orientation [oriDeg] in degrees.
%
% INPUT:
%     [nOri] : # of orientations bins 
%   [oriDeg] : desired orientation. Horizontal = 0, and vertical = 90. 
% [opt       : option structure
% .multiBin] : Use multiple bins when orientation falls on bin edge.
%              Default=0=use only single bin.
% OUTPUT:
% [oriIdx] : cell array of orientation index in coefficient cell arrays
%            that contains the orientation [oriDeg]. When the desired
%            orientation [oriDeg] is right on the boundary between two
%            orientation bins, it will give the index of both bins.
% [oriRng] : range of orientations that the bin index cover.
%
% SEE ALSO: crvOri, crvLet, crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Check Options
%--------------------
optDef.multiBin=0;
optRng.multiBin=[0,1];
if nargin<3
  opt=defaultOpt(optDef,optDef,optRng);
else
  opt=defaultOpt(opt,optDef,optRng);
end  



% Circularize the orientation in deg
%--------------------
oriDeg=mod(oriDeg,360);
nOriIn=length(oriDeg);



% Generate orientation bins
%--------------------
oriBinEdge=mod(135-linspace(0,360,nOri+1),360)';
oriBin=[oriBinEdge(2:end),oriBinEdge(1:end-1)];
zIdx=find(oriBinEdge==0);
oriBin(zIdx,2)=360;



% Determine which bin the data fall in
%--------------------
oriIdx=cell(nOriIn,1);
zinIdx=find(oriDeg==0);
nziIdx=find(oriDeg~=0);

% for data that are NOT 0=360 deg
if opt.multiBin
  for ii=nziIdx
    binIdx=oriDeg(ii)>=oriBin(:,1) & oriDeg(ii)<=oriBin(:,2);
    oriIdx{ii}=find(binIdx)';
  end  % for ii
else
  for ii=nziIdx  % for data that are NOT 0=360 deg
    binIdx=oriDeg(ii)>=oriBin(:,1) & oriDeg(ii)<oriBin(:,2);
    oriIdx{ii}=find(binIdx)';
  end  % for ii
end  % if opt.multiBin

% for data that ARE 0=360 deg
if opt.multiBin
  for ii=zinIdx
    Idx000=  0>=oriBin(:,1) &   0<=oriBin(:,2);
    Idx360=360>=oriBin(:,1) & 360<=oriBin(:,2);
    oriIdx{ii}=wrev(find(Idx000 | Idx360)');
  end
else
  for ii=zinIdx
    Idx000=  0>=oriBin(:,1) &   0<oriBin(:,2);
    Idx360=360>=oriBin(:,1) & 360<oriBin(:,2);
    oriIdx{ii}=wrev(find(Idx000 | Idx360)');
  end
end


 
% Generate orientation range
%--------------------
if nargout>1
  oriRng=nan(nOriIn,2);
  for ii=nziIdx  % for data that are NOT 0=360 deg
    oriRng(ii,:)=[oriBin(oriIdx{ii}(end),1),oriBin(oriIdx{ii}(1),2)];
  end
  for ii=zinIdx
    oriRng(ii,:)=[oriBin(oriIdx{ii}(end),1),oriBin(oriIdx{ii}(1),2)];
  end
end  % if nargout


