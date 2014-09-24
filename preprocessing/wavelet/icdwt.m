% icdwt.m
%
% Wrapper function for NGK's 1D dual-tree complex wavelet code.
% The wavelet set (near_sym_a, qshift_a) is hardwired.
% Usage : x = icdwt(w, L)
%
% Written by : Justin Romberg
% Created : 12/5/2000

function x = icdwt(w, L)

rw = 0;
if (size(w,1) == 1)
  rw = 1;
  w = w.';
end

Lx = log2(length(w));
% build the I matrix
I = zeros(L, 2);
for ll = Lx-L+1:Lx
  I(ll-Lx+L,:) = [2^ll 1];
end

C = zeros(2*length(w),1);
[cm,V] = icwtband(w(1:2^(Lx-L)), I, L, 'l');
C(V) = cm;
for ll = Lx:-1:(Lx-L+1)
  inds = 2^(ll-1)+1:2^ll;
  [cm,V] = icwtband(w(inds), I, Lx-ll+1, 'h');
  C(V) = cm;
end

x = dtwaverec(C, I, 'near_sym_a', 'qshift_a');
if (rw == 1)
  x = x';
end
