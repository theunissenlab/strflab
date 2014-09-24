function [cfs]=cwtZeroLev(cfs,lev,c)
%function [cfs]=cwtZeroLev(cfs,lev,c)
%
% Set to zero, or a constant [c], all the coefficients
% of the CWT at the specified scale level. The allow levels
% are determined by the # of levels use for to compute the 
% CWT. If N levels are used, there will be total of N+1
% levels, where the (N+1)th level is the low pass level.
% If not output arguments are specified, it will display
% the modified coefficiet subimages.
%
% INPUT:
%  [cfs] = cell array of coefficients at each scale level, with
%          orientated coeff of each level is stored as a structure
%          array inside the cell array of scales.
%  [lev] = Array of levels to zero out.
%    [c] = a constant, for which all CWT coefficients will
%          be set to. Default [c]=0.
% OUTPUT:
%  [cfs] = same as input [cfs] but with modified coeffs.
%
% SEE ALSO: cwtZeroOri, cwtZeroBlock, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input & Inits
%--------------------
if nargin<3
  c=0;
end
lev=unique(lev);
nLev=length(lev);
cfLev=length(cfs)-1;
if any(lev>cfLev+1)
  error('zLevCoef >< [lev] exceeds # of existing levels');
end


% Zero out coefficients lowpass filter level
%--------------------
if lev(end)>cfLev
  zblk=repmat(c,size(cfs{lev(end)}));
  cfs{lev(end)}=zblk;
  lev(end)=[];
end


% Zero out coefficients at [lev]
%--------------------
for ii=lev
  zblk=repmat(c,size(cfs{ii}.p75));
  cfs{ii}.p75=zblk;
  cfs{ii}.p45=zblk;
  cfs{ii}.p15=zblk;
  cfs{ii}.n15=zblk;
  cfs{ii}.n45=zblk;
  cfs{ii}.n75=zblk;
end  % for ii


% Show modified CWT coeff as subimages.
%--------------------
if nargout<1
  cwtShowStr(cfs);
end

