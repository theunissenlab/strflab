function [qcf]=cwtStr2qcf(cf)
%function [qcf]=cwtStr2qcf(cf)
%
% Convert the sub-image coeff structure [cf] back to the Dual
% Tree CWT coeff vector with size structure [qcf]. This is 
% needed to compute the inverse CWT.
%
% INPUT:
%   [cf] = CWT coeffs in sub-image structure. cell array of 
%          struct of coeffs at each scale level. Orientation 
%          of each level is stored as a structure. So that 
%            [cfStr{N}.p15] = Nth scale level at +15 deg.
%            [cfStr{M}.n45] = Mth scale level at -45 deg.
%          The highest level is the low pass level and has
%          no orientation structure.
%            [cfStr{end}] = low pass filter level. 
% OUTPUT:
%  [qcf] = CWT coeffs in dual tree coeff vector structure
%          containing [.coef] = complex wavelet coefficients 
%          and size of each level, and [.siz] size of subimages
%          at each level.
%
% SEE ALSO: cwtQcf2str, cwtDT, cwtInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Get data size & Init
%--------------------
lev=length(cf)-1;
blkDim=zeros(lev+1,2);
for ss=1:lev
  blkDim(ss,:)=size(cf{ss}.p75);
end
blkDim(lev+1,:)=size(cf{lev+1});
ncoef=2*sum((blkDim(1:end-1,1).*blkDim(1:end-1,2))*6)+blkDim(end,1).*blkDim(end,2);
qcf.coef=zeros(ncoef,1);
qcf.siz=2*flipud(blkDim(1:end-1,:));


% Fill in coefficients of lowpass filter level
%--------------------
[cIns,cIdx]=icwtband2(cf{lev+1},qcf.siz,lev,'l','real');
qcf.coef(cIdx)=cIns;


% Fill in coefficients of oriented bands
%--------------------
for ss=1:lev
  subI=cat(3,cf{ss}.p15,cf{ss}.n15, ...
             cf{ss}.p75,cf{ss}.n75, ...
             cf{ss}.p45,cf{ss}.n45);
  if isreal(subI)
    subI=complex(subI);
  end
  [cIns,cIdx]=icwtband6(subI,qcf.siz,ss);
  qcf.coef(cIdx)=cIns;
end  % for ii


