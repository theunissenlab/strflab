function [ccoef,imgSize,nLev]=crvLet(img,useWavlet,nLev,nAngle)
%function [ccoef,imgSize,nLev]=crvLet(img,useWavlet,nLev,nAngle)
%
% Takes curvelet transform of an image to get the curvelet coefficients.
%
% INPUT:
%       [img] : A matrix representing the luminance value at each pixel
%               of an image.
% [useWavlet] : 1 = USE wavelet at the finest scale = DEFAULT.
%               0 = DON'T use wavelet at the finest scale, use curvelet. 
%    [nLev] : # of scales for curvelet transform.
%               [DEFAULT = ceil(log2(min(imgSize(1),imgSize(2)))-3)].
%    [nAngle] : # of angles at the 2nd coarse scale. [DEFAULT = minimum 
%               = 8]. Must be multiples of 4.
% OUTPUT:
%     [ccoef] : cell array of curvelet coefficients at each scale [s],
%               orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
%   [imgSize] : size of image [img].
%    [nLev] : # of scales for curvelet transform.
%
% SEE ALSO: crvInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
imgSize=size(img); 
if nargin<4 | isempty(nAngle);
    nAngle=8;
end
if nargin<3 | isempty(nLev);
    nLev=ceil(log2(min(imgSize(1),imgSize(2)))-3);
end
if nargin<2 | isempty(useWavlet);
    useWavlet=1;  % Whether to use wavelet at the finest scale
end


% Taking Curvelet Transform
%--------------------
% 2nd argument 0 = complex wavelet
ccoef=fdct_wrapping(double(img),0,useWavlet+1,nLev,nAngle);

