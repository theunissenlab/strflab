function [cfMat,strSiz]=mov2cont(mov,opt)
%function [cfMat,strSiz]=mov2cont(mov,opt)
%
% Decompose each frame of a movie into contourlet decompositions.
% So each frame becomes a row vector of contourlet coeffs.
%
% INPUT:
%   [mov] : 3D matrix of size p x p x T. Where p=# of pixels for
%           each frame, and T=# of frames.
%   [opt] : option structure for contourlet transform.
% 
% OUTPUT:
%  [cfMat] : 2D matrix of ~3(p^2) x T, where each column of 
%            [mvccf] is a set of curvelet coeff.
% [strSiz] : Curvelet coeff structure for house keeping and
%            inversion of the features.
%
% SEE ALSO: 
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jul 2007)
%
%====================

% Get movie size/dimension
%--------------------
movSiz=size(mov);
nFrame=movSiz(3);
frmSiz=movSiz(1:2);


% Set Default Options
%--------------------
optDef.nLev=ceil(log2(min(frmSiz)))-1;
optRng.nLev=[1,16];
if nargin<2
  opt=optDef;
else
  opt=defaultOpt(opt,optDef,optRng);
end


% Compute curvelet transform for 1 frame
%--------------------
cfStr=contLet(mov(:,:,1),opt.nLev);
[cf1mat,strSiz]=pdfb2vec(cfStr);
matLen=length(cf1mat);

if opt.verbose, disp('Preprocessing...'); end
% Initialize output coeff matrix
%--------------------
cfMat=single(zeros(nFrame,matLen));
cfMat(1,:)=cf1mat;


% Compute curvelet transform for all frames
%--------------------
for ii=2:nFrame
  cfStr=contLet(mov(:,:,ii),opt.nLev, opt.pyr, opt.dir);
  cfMat(ii,:)=pdfb2vec(cfStr);
	if opt.verbose
		if mod(ii, 50) == 0, fprintf('.'), end
		if mod(ii, 1000) == 0, fprintf(' %d frames done.\n', ii), end
	end
end


