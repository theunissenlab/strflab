function [kScale]=cwtKeyScale(emap,keyPos)
%function [kScale]=cwtKeyScale(emap,keyPos)
%
% Given and accumulated energy map [emap] with its local
% optimal points as the location of the key points [keyPos],
% this function determines the scale of the key points.
%
% INPUT:
%   [emap] = 2D matrix representing an accumulated energy
%            map of an image obtained via cwtEnergy.m
% [keyPos] = Nx2 matrix, of key point locations, which are
%            the local maxima of the [emap].
% OUTPUT:
% [kScale] = Nx1 vector of key point scales.
%
% SEE ALSO: cwtEnergy, plotCircle
%
% By Michael Wu  --  waftingpetal@yahoo.com (May 2007)
%
% ====================


% Get Key Point Info
%--------------------
nKey=size(keyPos,1);
mapSiz=size(emap);
klen=zeros(1,4);
kScale=zeros(nKey,1);


for ii=1:nKey
  % Get Key Point Info
  %--------------------
  col.idx=keyPos(ii,1);
  row.idx=keyPos(ii,2);
 
  row.inc=row.idx:mapSiz(2);
  klen(1)=length(row.inc);
  row.dec=row.idx:-1:1;
  klen(2)=length(row.dec);
  row.con=repmat(row.idx,1,mapSiz(2));
  
  col.inc=col.idx:mapSiz(1);
  klen(3)=length(col.inc);
  col.dec=col.idx:-1:1;
  klen(4)=length(col.dec);
  col.con=repmat(col.idx,1,mapSiz(1));

  % Get Projection Direction
  %--------------------
  minLen=min([klen,ceil(mapSiz/4)]);
  cutIdx=zeros(minLen,8);
  
  cutIdx(:,1)=sub2ind(mapSiz,row.inc(1:minLen),col.con(1:minLen))';
  cutIdx(:,2)=sub2ind(mapSiz,row.dec(1:minLen),col.con(1:minLen))';
  cutIdx(:,3)=sub2ind(mapSiz,row.con(1:minLen),col.inc(1:minLen))';
  cutIdx(:,4)=sub2ind(mapSiz,row.con(1:minLen),col.dec(1:minLen))';
  cutIdx(:,5)=sub2ind(mapSiz,row.inc(1:minLen),col.inc(1:minLen))';
  cutIdx(:,6)=sub2ind(mapSiz,row.inc(1:minLen),col.dec(1:minLen))';
  cutIdx(:,7)=sub2ind(mapSiz,row.dec(1:minLen),col.inc(1:minLen))';
  cutIdx(:,8)=sub2ind(mapSiz,row.dec(1:minLen),col.dec(1:minLen))';
  
  % Form Gradient Projections
  %--------------------
  keyCut=emap(cutIdx);
  difCut=[zeros(1,8);diff(keyCut,1,1)];
  %difCut(difCut>0)=0;
  cutAxes=cumsum(ones(size(difCut)),1)-1;
  difCut(:,5:8)=difCut(:,5:8)/sqrt(2);
  cutAxes(:,5:8)=cutAxes(:,5:8)*sqrt(2);

  % Interpolate Gradient Projections
  %--------------------
  intPt=linspace(cutAxes(1,1),cutAxes(end,1),100);
  intVal=zeros(100,8);
  for jj=1:8
    intVal(:,jj)=interp1(cutAxes(:,jj),difCut(:,jj),intPt,'pchip');
  end
  
  % Extract Key Point Scale
  %--------------------
  [dummy,minIdx]=min(mean(intVal'));
  kScale(ii)=intPt(min(minIdx));
  
end  % for ii


