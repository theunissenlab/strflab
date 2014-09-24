function [cfmat,strSiz]=cwtStr2mat(cfStr)
%function [cfmat,strSiz]=cwtStr2mat(cfStr)
%
% Convert a subband image [cfStr] of the coeff of CWT into a 
% vector of coefficients [cfmat].
%
% INPUT:
%  [cfStr] = cell array of each scale, containing struct arrays 
%            of CWT coeff at different orientation subband. 
%            Orientation of each level is stored as a structure.
%            So that 
%              [cfStr{N}.p15] = Nth scale level at +15 deg.
%              [cfStr{M}.n45] = Mth scale level at -45 deg.
%            The highest level is the low pass level and has
%            no orientation structure.
%              [cfStr{end}] = low pass filter level. 
% OUTPUT:
%  [cfmat] = complex coeffs in a row vector.
% [strSiz] = coeff structure size and index in [cfmat] with 
%            fields [.siz] that shows the size of the original
%            block of coefficients. And field [.idx] contains
%            where in [cfmat] that block of coeff is stored.
% 
% SEE ALSO: cwtMat2str, cwtQcf2str, cwtDT, cwtInv, cwtShowStr
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Check input and Inits
%--------------------
lev=length(cfStr)-1;
strSiz=cell(1,lev+1);
blkSiz=zeros(1,lev);


% Get data size & Init [cfmat]
%--------------------
for ss=1:lev
  strSiz{ss}.siz=size(cfStr{ss}.p75);
  blkSiz(ss)=prod(strSiz{ss}.siz);
end
fineSiz=size(cfStr{lev+1});
blkSiz=repmat(blkSiz,6,1);
blkVec=[blkSiz(:);prod(fineSiz)];
cfIdx2=cumsum(blkVec(:));
cfIdx1=[0;cfIdx2(1:end-1)]+1;
cfIdx=[cfIdx1,cfIdx2];
cfmat=zeros(1,cfIdx(end,end));


% Convert cwt directional subband structure to a vector of coeff
%--------------------
for ss=1:lev
  strSiz{ss}.idx=cfIdx(1:6,:);
  cfIdx(1:6,:)=[];
  insIdx(1)=strSiz{ss}.idx(1,1);
  insIdx(2)=strSiz{ss}.idx(6,2);
  insLen=insIdx(2)-insIdx(1)+1;
  cfmat(insIdx(1):insIdx(2))=reshape(cell2mat(struct2cell(cfStr{ss})'),1,insLen);
end


% Convert cwt Lo-Pass band structure to a vector of coeff
%--------------------
strSiz{lev+1}.siz=fineSiz;
strSiz{lev+1}.idx=cfIdx;
cfmat(cfIdx(1):cfIdx(2))=reshape(cfStr{lev+1},1,prod(fineSiz));


