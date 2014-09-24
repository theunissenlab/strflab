function [strf,w]=lin1_2OrdPak(strf)
%function [strf, w] = linPak(strf)
% Combines weights and biases of a linear model into one weight vector.
%	
% Takes a strf data structure [strf] and  combines them
% into a single row vector [w]
%
% INPUTS:
%	
%	strf = a linear model strf structure (see linInit)
%
%
% OUTPUTS:
% strf = unmodified strf structure
%	 w = parameter vector for [strf] (ex. w=strfPak(strf))
%
%
%
%	SEE ALSO
%   linInit, linUnpak
%
%(Some code modified from NETLAB)

w=[strf.w1(:)' strf.w2(:)'];
if isfield(strf,'b1')
  w=[w,strf.b1];
end

