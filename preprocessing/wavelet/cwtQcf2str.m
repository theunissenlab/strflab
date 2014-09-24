function [cfStr]=cwtQcf2str(qcf)
%function [cfStr]=cwtQcf2str(qcf)
%
% Extract the sub-images structure of wavelet coefficients at
% each scale level and store in [cfStr]. If no output is 
% requested, it will display the coefficients in the subimages.
%
% INPUT:
%   [qcf] = structure containing [.coef] = complex wavelet
%           coefficients and size of each level, and [.siz]
%           size of subimages at each level.
% OUTPUT:
% [cfStr] = cell array of struct of coeffs at each scale level.
%           Orientation of each level is stored as a structure.
%           So that 
%             [cfStr{N}.p15] = Nth scale level at +15 deg.
%             [cfStr{M}.n45] = Mth scale level at -45 deg.
%           The highest level is the low pass level and has
%           no orientation structure.
%             [cfStr{end}] = low pass filter level. 
%
% SEE ALSO: cwtDT, cwtInv, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check Input & Inits
%--------------------
lev=size(qcf.siz,1);
cfStr=cell(1,lev+1);


% Extract oriented Scale levels
%--------------------
for ss=1:lev
  band=cwtband6(qcf.coef,qcf.siz,ss);
  cfStr{ss}.p75=band(:,:,3);
  cfStr{ss}.p45=band(:,:,5);
  cfStr{ss}.p15=band(:,:,1);
  cfStr{ss}.n15=band(:,:,2);
  cfStr{ss}.n45=band(:,:,6);
  cfStr{ss}.n75=band(:,:,4);
end  % for ii


% Extract low pass level
%--------------------
cfStr{lev+1}=cwtband2(qcf.coef,qcf.siz,lev,'l','real');


% Show CWT coeff as subimages.
%--------------------
if nargout<1
  cwtShowStr(cfStr);
end

