function [simdex,strSim,nCoef]=cwtssim(cwt1,cwt2,lev)
%function [simdex,strSim,nCoef]=cwtssim(cwt1,cwt2,lev)
%
% Structural Similarity base on Dual Tree CWT coeffs.
%
% INPUT:
% [cwt1] = CWT coefficient structures. Field [.coef] contains
% [cwt2]   he real coefficients, and fiend [.siz] contains 
%          size of sub-image at each level.
%  [lev] = The scale level to compare the similarity.
% OUTPUT:
% [strSim] = The similarity measure of each level compared.
% [simdex] = average sMean of [strSim].
%  [nCoef] = # of coefs at each of [lev] used to weight the
%            correlation coefficients.
%
% SEE ALSO: pairssim, dispPairCmp
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Initialization
%--------------------
nLev=length(lev);
%cwt1.coef=full(cwt1.coef);
%cwt2.coef=full(cwt2.coef);
si1=cwtCoef2SubIm(cwt1);
si2=cwtCoef2SubIm(cwt2);
strSim=zeros(nLev,1);
nCoef=strSim;


% Compute structural similarity at each scale.
%--------------------
for ll=1:nLev
  if isstruct(si1{lev(ll)})
  % if level has oriented coeff
    c1=cell2mat(struct2cell(si1{lev(ll)}));
    c2=cell2mat(struct2cell(si2{lev(ll)}));
    nCoef(ll)=length(c1(:));
    
    cpCorr=c1.*conj(c2);
    cpVar=abs(c1).^2+abs(c2).^2;
    tmp=2*abs(cpCorr)./cpVar;
    strSim(ll)=mean(tmp(:));
  else
  % if level has NO oriented coeff
    c1=si1{lev(ll)}(:);
    c2=si2{lev(ll)}(:);
    nCoef(ll)=length(c1);
    
    tmp=abs(corrcoef(c1,c2));
    strSim(ll)=tmp(1,2);
  end  % if isstruct    
end  % for ll


% Compute Weighted average of structural similarity
%--------------------
simdex=mean(strSim);
%simdex=sum(strSim./nCoef)./sum(1./nCoef);



