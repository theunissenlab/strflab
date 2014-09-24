function [mvccf,strSiz]=mov2crv(mov,opt)
%function [mvccf,strSiz]=mov2crv(mov,opt)
%
% Decompose each frame of a movie into curvelet decompositions.
% So each frame becomes a row vector of curvelet coeffs.
%
% INPUT:
%   [mov] : 3D matrix of size p x p x T. Where p=# of pixels for
%           each frame, and T=# of frames.
%   [opt] : option structure for curvelet transform.
% [.]
% OUTPUT:
%  [mvccf] : 2D matrix of ~3(p^2) x T, where each column of 
%            [mvccf] is a set of curvelet coeff.
% [strSiz] : Curvelet coeff structure for house keeping and
%            inversion of the features.
%
% SEE ALSO: crvMat2str, crvStr2mat, crvLet, crvInv.
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Get movie size/dimension
%--------------------
movSiz=size(mov);
nFrame=movSiz(3);
frmSiz=movSiz(1:2);


% Set Default Options
%--------------------
optDef.nLev=ceil(log2(min(frmSiz(1),frmSiz(2)))-3);
optDef.nAngle=8;
optDef.useWavelet=0;
optRng.nLev=[1,16];
optRng.nAngle=[8:4:32];
optRng.useWavelet=[0,1];
if nargin<2
  opt=optDef;
else
  opt=defaultOpt(opt,optDef,optRng);
end


% Compute curvelet transform for 1 frame
%--------------------
ccf=crvLet(mov(:,:,1),opt.useWavelet,opt.nLev,opt.nAngle);
[ccm,strSiz]=crvStr2mat(ccf,frmSiz);
cfLen=length(ccm);

if opt.verbose, disp('Preprocessing...'); end
% Initialize output coeff matrix
%--------------------
mvccf=single(zeros(nFrame,cfLen));
mvccf(1,:)=ccm;


% Compute curvelet transform for all frames
%--------------------
for ii=2:nFrame
  ccf=crvLet(mov(:,:,ii),opt.useWavelet,opt.nLev,opt.nAngle);
  mvccf(ii,:)=crvStr2mat(ccf,frmSiz);
	if opt.verbose
		if mod(ii, 50) == 0, fprintf('.'), end
		if mod(ii, 1000) == 0, fprintf(' %d frames done.\n', ii), end
	end
end


