function [cf,filtName]=contLet(im,nLev,pyr,dir)
%function [cf,filtName]=contLet(im,nLev,pyr,dir)
%
% Compute the coeff of contourlet transform.
%
% INPUT:
%   [im] = 2-D array of pixel values of an image.
% [nLev] = # of scale level of contourlet transform
%  [pyr] = scalar in [1,9] specifying the pyramidal filter type 
%              1='9-7'      2='5-3'      3='Burt'
%              4='pkva'     5='haar'     6='db2'
%              7='coif1'    8='sym2'     9='dmey'=Default.
%  [dir] = scalar in [1,7] specifying the directional filter type
%              1='haar'     2='5-3'=Default     
%              3='9-7'      4='pkva'     5='pkva6'
%              6='pkva8'    7='pkva12'
% OUTPUT:
%   [cf] = cell array of contourlet coefficients at each scale [s],
%          orientation [w], and position [x,y] as cf{s}{w}[x,y].
%          cf{1}=low pass level with no orientation. Increasing 
%          [s] gives finer scale.
% [filtName.pyr] = Name of pyramidal filter
% [filtName.dir] = Name of directional filter
%
% SEE ALSO: contInv
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jun 2007)
%
% ====================


% Init Parameters
%--------------------
lev=[1 2 3 3 4 4 5 5 6 6 7 7];
if nargin<3; pyr=9; end
if nargin<4; dir=2; end 

switch pyr % Possible pyramidal filter names
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

switch dir % Possible directional filter names
  case 1; dirFilt='haar';
  case 2; dirFilt='5-3'; % McClellan transformed of 5-3 filters
  case 3; dirFilt='9-7'; % by Cohen and Daubechies
  case 4; dirFilt='pkva';
  case 5; dirFilt='pkva6';
  case 6; dirFilt='pkva8';
  case 7; dirFilt='pkva12';
  otherwise dirFilt='5-3';
end


% Contourlet transform
%--------------------
cf=pdfbdec(double(im),pyrFilt,dirFilt,lev(1:nLev));
filtName.pyr=pyrFilt;
filtName.dir=dirFilt;

