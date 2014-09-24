function strf = nnunpak(strf, w)
%NNUNPAK Separates weights vector into weight and bias matrices. 
%
%	Description
%	STRF = NNUNPAK(STRF, W) takes an network data structure STRF and  a
%	weight vector W, and returns a network data structure identical to
%	the input network, except that the first-layer weight matrix W1, the
%	first-layer bias vector B1, the second-layer weight matrix W2 and the
%	second-layer bias vector B2 have all been set to the corresponding
%	elements of W.
%
%Some code taken from Netlab

if strf.nWts ~= length(w)
  error('Invalid weight vector length')
end

nin = strf.nIn;
nhidden = strf.nHidden;

mark1 = nin*length(strf.delays);

for ii=1:nhidden
  strf.w1(:,ii,:) = reshape(w(1:mark1), [nin 1 length(strf.delays)]);
  w=w(1,mark1+1:end);
end

mark2 = nhidden;
strf.b1 = reshape(w(1,1:mark2), 1, nhidden);

mark3 = mark2 + nhidden;
strf.w2 = reshape(w(mark2 + 1: mark3), nhidden, 1);

mark4 = mark3 + 1;
strf.b2 = reshape(w(mark3 + 1: mark4), 1, 1);
