function [strf,w] = midPak(strf)
%MIDPAK	Combines weights and biases into one weights vector.
%
%	Description
%	W = MIDPAK(STRF) takes a strf data structure STRF and  combines them
%	into a single row vector W.
%
%
%	Usage:
%	
%	W = MIDPAK(STRF)
%
%
%	Inputs:
%	
%	STRF,	A STRF structure
%
%
%	Outputs:
%
%	W,	Parameter vector for STRF, (ex. w=strfPak(strf))
%
%
%
%	See also
%	MID, MIDUNPAK, MIDFWD, MIDERR, MIDGRAD
%
%(Some code modified from NETLAB)

% errstring = consist(net, 'mid');
% if ~errstring
%   error(errstring);
% end

w=strf.w1(:)';
if isfield(strf,'b1')
  w=[w,strf.b1];
end
