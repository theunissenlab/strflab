function [strf,w]=strfPak(strf)
%function [strf,w]=strfPak(strf)
%
% Combines weights and biases into one weights vector.
%
% The facility to switch between these two representations for the network
% parameters is useful, for example, in training a strf by error
% function minimization, since a single vector of parameters can be
% handled by general-purpose optimization routines.
%
% INPUT:	
% [strf] = strf model structure
%
% OUTPUT:
% [strf] = strf model structure
%    [w] = Parameter vector for STRF, (ex. w=strfPak(strf))
%
%
%(Some code modified from NETLAB)

pakstr=[strf(1).type,'Pak'];
w=zeros(length(strf),strf(1).nWts);

for ii=1:length(strf)
  [strf(ii),w(ii,:)]=feval(pakstr,strf(ii));
end
