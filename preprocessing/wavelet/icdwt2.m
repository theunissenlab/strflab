% icdwt2.m
%
% Wrapper function for NGK's 2D dual-tree complex wavelet code
% Inverse 2D transform (synthesis)
% The wavelet set (near_sym_a, qshift_a) is hardwired in right now.
% Usage : x = icdwt2(w1, w2, L)
%
% Written by : Justin Romberg
% Created : 1/30/2001

function x = icdwt2(w1, w2, L)

N = size(w1,1);
Lx = log2(N);
% build the S matrix
S = zeros(L, 2);
for ll = Lx-L+1:Lx
  S(ll-Lx+L,:) = [2^ll 2^ll];
end

C = zeros(4*N*N,1);
[cm, V] = icwtband2([w1(1:2^(Lx-L),1:2^(Lx-L)); ...
  w2(1:2^(Lx-L),1:2^(Lx-L))], S, L, 'l');
C(V) = cm;
for ll = Lx:-1:(Lx-L+1)
  k = 2^(ll-1);
  [cm, V] = icwtband2([w1(1:k,k+1:2*k); w2(1:k,k+1:2*k)], S, Lx-ll+1, 'v');
  C(V) = cm;
  [cm, V] = icwtband2([w1(k+1:2*k,1:k); w2(k+1:2*k,1:k)], S, Lx-ll+1, 'h');
  C(V) = cm;
  [cm, V] = icwtband2([w1(k+1:2*k,k+1:2*k); w2(k+1:2*k,k+1:2*k)], ...
      S, Lx-ll+1, 'd');
  C(V) = cm;
end

x = dtwaverec2(C, S, 'near_sym_a', 'qshift_a');


    
