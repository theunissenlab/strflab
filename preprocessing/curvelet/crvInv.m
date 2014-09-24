function [img,nScale]=crvInv(ccoef,imSiz)
%function [img,nScale]=crvInv(ccoef,imSiz)
%
% Takes inverse curvelet transform of curvelet coefficients to recover an
% image.
%
% INPUT:
%  [ccoef] : cell array of curvelet coefficients at each scale [s],
%            orientation [w], and position [x,y]. ccoef{s}{w}[x,y].
%  [imSiz] : [mRow,nCol] = # of rows (height of the image [img], and
%            # of columns (width of the image [img].
%
% OUTPUT:
%    [img] : A [mRow] by [nCol] matrix representing the luminance value 
%            at each pixel of the reconstructed image.
% [nScale] : # of scales for curvelet transform.
%
% SEE ALSO: crvLet
%
% By Michael Wu  --  waftingpetal@yahoo.com (Sep 2006)
%
%====================


% Init Parameters
%--------------------
if nargin<2
  coefLen=length(ccoef{end});
  if coefLen~=1
    error('crvInv >< must provide [imgSiz]');
  else
    imSiz=size(ccoef{end}{end});
  end
end


% Inverse Curvelet Transform
%--------------------
nScale=length(ccoef);
img=ifdct_wrapping(ccoef,0,imSiz(1),imSiz(2));


% Display image reconstruction
%--------------------
if nargout<1
  hf=figure;
  imagesc(real(img));
end  % if nargout
