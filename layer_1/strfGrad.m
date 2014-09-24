function [strf,g]=strfGrad(strf,datIdx)
%function [strf,g]=strfGrad(strf,datIdx)
%
% Evaluates strf error gradient for generic optimizers
%
% INPUT:	
%   [strf] = strf model structure
% [datIdx] = a vector containing indices of the samples to be used in the 
%            fitting.
%
% OUTPUT:
%   [strf] = strf model structure
%      [g] = value of the gradient of the STRF's error function for the
%            indexed data
%
%
%(Some code modified from NETLAB)


gradstr = [strf.type,'Grad'];
[strf,g]=feval(gradstr,strf,datIdx);
