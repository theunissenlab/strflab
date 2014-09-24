function [cfs]=cwtZeroOri(cfs,ori,c)
%function [cfs]=cwtZeroOri(cfs,ori,c)
%
% Set to zero, or a constant [c], all the coefficients
% of the CWT for the specified orientation [ori] at all level. 
% Ori is specified by:
%      'v' = vertical orientation   = (+/-) 75 degrees
%      'h' = horizontal orientation = (+/-) 15 degrees
%      'd' = diagonal orientation   = (+/-) 45 degrees
%      'l' = low pass level (no orientation).
% If not output arguments are specified, it will display
% the modified coefficiet subimages.
%
% INPUT:
%  [cfs] = cell array of coefficients at each scale level, with
%          orientated coeff of each level is stored as a structure
%          array inside the cell array of scales.
%  [ori] = char array of orientations to zero out.
%    [c] = a constant, for which all CWT coefficients will
%          be set to. Default [c]=0.
% OUTPUT:
%  [cfs] = same as input [cfs] but with modified coeffs.
%
% SEE ALSO: cwtZeroLev, cwtZeroBlock, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input
%--------------------
if nargin<3
  c=0;
end
ori=unique(lower(ori));
cfLev=length(cfs)-1;
if ~isempty(setdiff(ori,'vhdl'))
  error('zroLevOri >< [ori] must be char array of ''v'', ''h'', ''d'' &/or ''l'' ');
end


% Zero out coefficients lowpass level
%--------------------
if ismember('l',ori)
  zblk=repmat(c,size(cfs{cfLev+1}));
  cfs{cfLev+1}=zblk;
  ori=setdiff(ori,'l');
end


% Zero out oriented coefficients
%--------------------
for ii=1:cfLev
  zblk=repmat(c,size(cfs{ii}.p75));
  if ismember('h',ori)
    cfs{ii}.p15=zblk;
    cfs{ii}.n15=zblk;
  end
  if ismember('v',ori)
    cfs{ii}.p75=zblk;
    cfs{ii}.n75=zblk;
  end
  if ismember('d',ori)
    cfs{ii}.p45=zblk;
    cfs{ii}.n45=zblk;
  end
end


% Show modified CWT coeff as subimages.
%--------------------
if nargout<1
  cwtShowStr(cfs);
end

