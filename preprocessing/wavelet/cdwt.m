% cdwt.m
%
% Wrapper function for NGK's 1D dual-tree complex wavelet code.
% The wavelet set (near_sym_a, qshift_a) is hardwired in right now.
% Usage : w = cdwt(x, L)
%
% Written by : Justin Romberg
% Created : 12/5/2000

function w = cdwt(x, L)

% make x a column vector
rw = 0;
if (size(x,1) == 1)
  x = x';
  rw = 1;
end

Lx = log2(length(x));

% this is what we usually use
[W,I] = dtwavedec(x, L, 'near_sym_a', 'qshift_a');
%[W,I] = dtwavedec(x, L, 'antonini', 'qshift_b');
w = zeros(size(x));
for ll = Lx:-1:(Lx-L+1)
  inds = 2^(ll-1)+1:2^ll;
  w(inds) = cwtband(W, I, Lx-ll+1, 'h');
end
w(1:2^(Lx-L)) = cwtband(W, I, L, 'l');

% if x was originally a row, return a row
if (rw == 1)
  w = w.';
end


