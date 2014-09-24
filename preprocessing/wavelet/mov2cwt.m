function [cfMat,strSiz]=mov2cwt(mov,opt)
%function [cfMat,strSiz]=mov2cwt(mov,opt)
%
% Decompose each frame of a movie into complex wavelets.
% So each frame becomes a row vector of wavelet coeffs.
%
% INPUT:
%   [mov] : 3D matrix of size p x p x T. Where p=# of pixels for
%           each frame, and T=# of frames.
%   [opt] : option structure for complex wavelet transform.
% [.]
% OUTPUT:
%  [cfMat] : 2D matrix of ~3(p^2) x T, where each column of [mvccf]
%            is a set of complex wavelet coeff.
% [strSiz] : Curvelet coeff structure for house keeping and inversion
%            of the features.
%
% SEE ALSO: 
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jun 2007)
%
%====================


% Set Default Options
%--------------------
optDef.nLev=3;
optRng.nLev=[1,16];

if nargin<1
  cfMat=optDef;
  return;
end

if nargin<2
  opt=optDef;
else
  opt=defaultOpt(opt,optDef,optRng);
end


% Get movie size/dimension
%--------------------
movSiz=size(mov);
nFrame=movSiz(3);
frmSiz=movSiz(1:2);


% Compute curvelet transform for 1 frame
%--------------------
cfStr=cwtDT(mov(:,:,1),opt.nLev);
[cf1mat,strSiz]=cwtStr2mat(cfStr);
matLen=length(cf1mat);

if opt.verbose, disp('Preprocessing...'); end
% Initialize output coeff matrix
%--------------------
cfMat=single(zeros(nFrame,matLen));
cfMat(1,:)=cf1mat;


% Compute curvelet transform for all frames
%--------------------
for ii=2:nFrame
  cfStr=cwtDT(mov(:,:,ii),opt.nLev);
  cfMat(ii,:)=cwtStr2mat(cfStr);
	if opt.verbose
		if mod(ii, 50) == 0, fprintf('.'), end
		if mod(ii, 1000) == 0, fprintf(' %d frames done.\n', ii), end
	end
end


