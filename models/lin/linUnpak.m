function strf=linUnpak(strf,w)
%function [strf] = linUnpak(strf, w)
% Separates weights vector into weight and bias matrices. 
%
% Takes a lin strf data structure [strf] and a
% weight array [w], and returns a lin data structure with the fields strf.w1, and
% strf.b1 set to the elements in w.  The bias term "b1" should be the LAST element
% of the array [w].
%
% INPUTS:	
%	[strf] = a linear model strf structure
%	   [w] = parameter vector for STRF, (ex. w=strfPak(strf))
%
%
% OUTPUTS:
%
%	[strf] = a strf structure with the fields w
%
%
%	See also
%	linInit, linPak, linFwd, linErr, linGrad
%
%
%(Some code modified from NETLAB)

if strf.nWts~=length(w)
  error('Invalid weight vector length')
end

nIn=strf.nIn;
strf.w1=reshape(w(1:nIn*length(strf.delays)),nIn,1,length(strf.delays));

if isfield(strf,'b1')
  strf.b1=reshape(w(nIn*length(strf.delays)+1),1,1);
end

strf.internal.compFwd = 1;
