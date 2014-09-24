function [cfStr]=cwtMat2str(cfmat,strSiz)
%function [cfStr]=cwtMat2str(cfmat,strSiz)
%
% Convert a row vector of CWT coeff [cfmat] back to a subimage
% structure [cfStr]. Essentially the inverse of cwtStr2mat.
% 
% INPUT:
%  [cfmat] = complex coeffs in a row vector.
% [strSiz] = coeff structure size and index in [cfmat] with 
%            fields [.siz] that shows the size of the original
%            block of coefficients. And field [.idx] contains
%            where in [cfmat] that block of coeff is stored.
% OUTPUT:
%  [cfStr] = cell array of each scale, containing struct arrays 
%            of CWT coeff at different orientation subband. 
%            Orientation of each level is stored as a structure.
%            So that 
%              [cfStr{N}.p15] = Nth scale level at +15 deg.
%              [cfStr{M}.n45] = Mth scale level at -45 deg.
%            The highest level is the low pass level and has
%            no orientation structure.
%              [cfStr{end}] = low pass filter level. 
%
% SEE ALSO: cwtStr2mat, cwtQcf2str, cwtDT, cwtInv, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check input and Inits
%--------------------
lev=length(strSiz)-1;
cfStr=cell(1,lev+1);


% Convert cwt row vector coeff back to directional subband structure
%--------------------
for ss=1:lev
  cfStr{ss}.p75=reshape(cfmat(strSiz{ss}.idx(1,1):strSiz{ss}.idx(1,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
  cfStr{ss}.p45=reshape(cfmat(strSiz{ss}.idx(2,1):strSiz{ss}.idx(2,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
  cfStr{ss}.p15=reshape(cfmat(strSiz{ss}.idx(3,1):strSiz{ss}.idx(3,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
  cfStr{ss}.n15=reshape(cfmat(strSiz{ss}.idx(4,1):strSiz{ss}.idx(4,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
  cfStr{ss}.n45=reshape(cfmat(strSiz{ss}.idx(5,1):strSiz{ss}.idx(5,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
  cfStr{ss}.n75=reshape(cfmat(strSiz{ss}.idx(6,1):strSiz{ss}.idx(6,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));
end


% Convert cwt row vector coeff back to the lo-pass band
%--------------------
ss=lev+1;
cfStr{ss}=reshape(cfmat(strSiz{ss}.idx(1,1):strSiz{ss}.idx(1,2)),strSiz{ss}.siz(1),strSiz{ss}.siz(2));


