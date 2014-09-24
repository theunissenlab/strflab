function [ccf]=crvMat2str(ccm,strSiz)
%function [ccf]=crvMat2str(ccm,strSiz)
%
% Converts a curvelet coeff 1D matrix back to a structure
% for ease of inversion.
%
% INPUT:
%      [ccm] : curvelet coeff matrix
%   [strSiz] : size of original [ccf] structure and where
%              the scale & orientation coeffs are located
%              in [ccm]. In strSiz{ss}{oo} is a structure
%              the [.siz] field give the origianl size of
%              the coeffs. And [.idx] field gives the index
%              location in [ccm] where the coeff are found.
% OUTPUT:
%      [ccf] : curvelet coeff structure.
%
% SEE ALSO: crvMat2str, crvLet, crvInv, mov2crv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Get curvelet coeff structure size info
%--------------------
nLev=length(strSiz);
ccf=cell(1,nLev);


% Get curvelet coeff structure size info
%--------------------
for ss=1:nLev
  nOri=length(strSiz{ss});
  ccf{ss}=cell(1,nOri);
  for oo=1:nOri
    pSiz=strSiz{ss}{oo}.siz;
    pIdx=strSiz{ss}{oo}.idx;
    ccf{ss}{oo}=reshape(ccm(pIdx(1):pIdx(2)),pSiz(1),pSiz(2));  
  end  % for oo
end  % for ss

