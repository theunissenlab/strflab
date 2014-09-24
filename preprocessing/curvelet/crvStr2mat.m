function [ccm,strSiz,overComp]=crvStr2mat(ccf,imSiz)
%function [ccm,strSiz,overComp]=crvStr2mat(ccf,imSiz)
%
% Converts a curvelet coeff structure to a 1D matrix.
%
% INPUT:
%      [ccf] : curvelet coeff structure.
%    [imSiz] : dimension of original image.
% OUTPUT:
%      [ccm] : curvelet coeff matrix
%   [strSiz] : size of original [ccf] structure and where
%              the scale & orientation coeffs are located
%              in [ccm]. In strSiz{ss}{oo} is a structure
%              the [.siz] field give the origianl size of
%              the coeffs. And [.idx] field gives the index
%              location in [ccm] where the coeff are found.
% [overComp] : level of overcompleteness = ratio of # of 
%              curvelet coeff to # of pixels.
%
% SEE ALSO: crvMat2str, crvLet, crvInv, mov2crv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Check input
%--------------------
nLev=length(ccf);
finestLen=length(ccf{nLev});
if finestLen==1
  useWavelet=true;
  imSiz=size(ccf{nLev}{finestLen});
else
  useWavelet=false;
end


% Initialization
%--------------------
strSiz=cell(1,nLev);
matIdx=0;
ccm=zeros(1,3*prod(imSiz));


% Convert curvelet coeff structure to matrix
%--------------------
for ss=1:nLev
  nOri=length(ccf{ss});
  strSiz{ss}=cell(1,nOri);
  for oo=1:nOri
    pSiz=size(ccf{ss}{oo});
    ppSiz=prod(pSiz);
    strSiz{ss}{oo}.siz=pSiz;
    strSiz{ss}{oo}.idx=[matIdx+1,matIdx+ppSiz];
    ccm(matIdx+1:matIdx+ppSiz)=reshape(ccf{ss}{oo},1,ppSiz);
    matIdx=matIdx+ppSiz;
  end  % for oo
end  % for ss
ccm(matIdx+1:end)=[];


% Compute overcompleteness
%--------------------
if nargout>2
  overComp=matIdx/prod(imSiz);
end

