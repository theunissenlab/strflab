function strf = midUnpak(strf, w)
%MIDUNPAK Separates weights vector into weight and bias matrices. 
%
%	Description
%	STRF = MIDUNPAK(STRF, W) takes a mid strf data structure STRF and  a
%	weight vector W, and returns a strfwork data structure identical to
%	the input strfwork, except that the first-layer weight matrix W1 and
%	the first-layer bias vector B1 have been set to the corresponding
%	elements of W.
%
%
%	Usage:
%	
%	W = MIDUNPAK(STRF)
%
%
%	Inputs:
%	
%	STRF,	A STRF structure
%	W,	Parameter vector for STRF, (ex. w=strfPak(strf))
%
%
%	Outputs:
%
%	STRF,	A STRF structure
%
%
%	See also
%	MID, MIDPAK, MIDFWD, MIDERR, MIDGRAD
%
%
%(Some code modified from NETLAB)

% Check arguments for consistency
% errstring = consist(strf, 'mid');
% if ~errstring
%   error(errstring);
% end

if strf.nWts~=length(w)
  error('Invalid weight vector length')
end

nIn=strf.nIn;
strf.w1=reshape(w(1:nIn*length(strf.delays)),nIn,1,length(strf.delays));

if isfield(strf,'b1')
  strf.b1=reshape(w(nIn*length(strf.delays)+1),1,1);
end

strf.internal.compFwd = 1;
