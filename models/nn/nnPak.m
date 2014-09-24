function [strf,w]= nnPak(strf)
%NNPAK	Converts neural network strf structure to vector of free parameters.
%
%	Description
%	W = NNPAK(STRF) takes a network data structure STRF and combines the
%	component weight matrices bias vectors into a single row vector W.
%	The facility to switch between these two representations for the
%	network parameters is useful, for example, in training a network by
%	error function minimization, since a single vector of parameters can
%	be handled by general-purpose optimization routines.
%
%	The ordering of the paramters in W is defined by
%	  w = [strf.w1(:)', strf.b1, strf.w2(:)', strf.b2];
%	 where W1 is the first-layer weight matrix, B1 is the first-layer
%	bias vector, W2 is the second-layer weight matrix, and B2 is the
%	second-layer bias vector.
%
%

%	Some code taken from Netlab


w=[];
for ii=1:strf.nHidden
  w = [w, reshape(strf.w1(:,ii,:),[1 strf.nIn*length(strf.delays)])];
end
w = [w, strf.b1];
w = [w, strf.w2(:)'];
w = [w, strf.b2];


