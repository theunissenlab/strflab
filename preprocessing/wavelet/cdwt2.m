% cdwt2.m
%
% Wrapper function for NGK's 2D dual-tree complex wavelet code
% Forward 2D transform (analysis)
% The wavelet set (near_sym_a, qshift_a) is hardwired in right now.
% Usage : [w1, w2] = cdwt2(x, L)
% w1 - subbands with directions
%       -----------------
%       |       |       |
%       |   X   |  +75  |
%       |       |       |
%       -----------------
%       |       |       |
%       |  +15  |  +45  |
%       |       |       |
%       -----------------
% w2 - subbands with directions
%       -----------------
%       |       |       |
%       |   X   |  -75  |
%       |       |       |
%       -----------------
%       |       |       |
%       |  -15  |  -45  |
%       |       |       |
%       -----------------
%
% Written by : Justin Romberg
% Created : 1/30/2001

function [w1, w2] = cdwt2(x, L)

Lx = log2(size(x,1));

[C,S] = dtwavedec2(x, L, 'near_sym_a', 'qshift_a');
w1 = zeros(size(x));
w2 = zeros(size(x));
for ll = Lx:-1:(Lx-L+1)
  k = 2^(ll-1);
  b75 = cwtband2(C, S, Lx-ll+1, 'v');
  w1(1:k,k+1:2*k) = b75(1:k,:);
  w2(1:k,k+1:2*k) = b75(k+1:2*k,:);
  b15 = cwtband2(C, S, Lx-ll+1, 'h');
  w1(k+1:2*k,1:k) = b15(1:k,:);
  w2(k+1:2*k,1:k) = b15(k+1:2*k,:);
  b45 = cwtband2(C, S, Lx-ll+1, 'd');
  w1(k+1:2*k,k+1:2*k) = b45(1:k,:);
  w2(k+1:2*k,k+1:2*k) = b45(k+1:2*k,:);
end
bl = cwtband2(C, S, L, 'l');
w1(1:2^(Lx-L),1:2^(Lx-L)) = bl(1:2^(Lx-L),:);
w2(1:2^(Lx-L),1:2^(Lx-L)) = bl(2^(Lx-L)+1:2^(Lx-L+1),:);

