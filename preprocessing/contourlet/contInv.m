function [im]=contInv(cf,fname)
%function [im]=contInv(cf,filtName)
%
% Compute the inverse contourlet transform from coeff to get the
% image back.
%
% INPUT:
%        [cf] = cell array of contourlet coefficients at each scale
%               [s], orientation [w], and position [x,y] as 
%               cf{s}{w}[x,y]. cf{1}=low pass level with no 
%               orientation. Increasing [s] gives finer scale.
% [fname.pyr] = Name of pyramidal filter.   Default='dmey'
% [fname.dir] = Name of directional filter. Default='5-3'
% OUTPUT:
%        [im] = 2-D array of pixel values of an image.
%
% SEE ALSO: contLet
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jun 2007)
%
% ====================


% Check input
%--------------------
if nargin<2
  pyrFilt='dmey';
  dirFilt='5-3';
else
%   pyrFilt=fname.pyr;
%   dirFilt=fname.dir;
% end

switch fname.pyr % Possible pyramidal filter names
  case 1; pyrFilt='9-7';
  case 2; pyrFilt='5-3';
  case 3; pyrFilt='Burt';
  case 4; pyrFilt='pkva'; % filters from the ladder structure
  % Standard Wavelets
  case 5; pyrFilt='haar';
  case 6; pyrFilt='db2';
  case 7; pyrFilt='coif1';
  case 8; pyrFilt='sym2';
  case 9; pyrFilt='dmey';
  otherwise pyrFilt='dmey';
end

switch fname.dir % Possible directional filter names
  case 1; dirFilt='haar';
  case 2; dirFilt='5-3'; % McClellan transformed of 5-3 filters
  case 3; dirFilt='9-7'; % by Cohen and Daubechies
  case 4; dirFilt='pkva';
  case 5; dirFilt='pkva6';
  case 6; dirFilt='pkva8';
  case 7; dirFilt='pkva12';
  otherwise dirFilt='5-3';
end

end
% Inverse Contourlet
%--------------------
im=pdfbrec(cf,pyrFilt,dirFilt);

