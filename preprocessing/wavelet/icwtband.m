function [cm,V] = icwtband(ym,L,level,passband,real_or_cplx);

% Function to determine where to insert a modified subband ym in C
% and to convert from complex to real format if ym is complex.
% It is necessary to use 'C(V) = cm;' to then do the insertion.
% (For large arrays C, this is much more efficient than copying
% the whole of C into and out of the function.)
%
% [cm, V] = icwtband(ym,L,level,passband,real_or_cplx);
%
%     ym -> The modified (real or complex) subband vector to be inserted.
%     L -> A bookkeeping matrix containing the size of the appropriate cwtbands
%     level -> The level at which ym is to be inserted
%     passband ->  'l' => Low-pass  
%                  'h' => Hi-pass
%     real_or_cplx ->  'real' => If ym is a real vector
%                      'cplx' => If ym is a complex vector. This is the DEFAULT.
%
%     cm -> The modified (real) subband vector for insertion into C
%     V -> The start-to-finish vector of points where cm should be inserted in C
%
% For example:  [cm,V] = icwtband(l,L,2,'l','real');
% or:           [cm,V] = icwtband(h3,L,3,'l','cplx');
%               C(V) = cm;  % to insert the modified vector into C.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000

%set the default for real_or_cplx
if nargin < 5, 
   real_or_cplx = 'cplx';
end 

[a b] = size(L);

%convert to REAL if needs be
if real_or_cplx == 'cplx';
   ym = c2q1d(ym);
%   warning('Now the output vector, cm, is a REAL vector, and should be inserted in C in this form.');
   cm = ym(:);
elseif real_or_cplx == 'real';
   cm = ym(:);
else
   error('Please enter the correct term for the appropriate cwtband. see help ICWTBAND');
end


if passband == 'l'; 
   if level ~= a;
      warning(sprintf('Only the Lo cwtband of the highest level is to be inserted, in this case level %d', a));
   end
   if all(size(ym) == L(1,:));		%check to see if the sizes match
      start = 1;
      finish = L(1,1)*L(1,2);		%the Lo band is always the first entry in L.
      V = [start:finish];
   else
      error('The size of ym and the selected passband and level do not agree.');
   end      
   
elseif passband == 'h';   
   size_seg = L(a+1-level,1)*L(a+1-level,2);   %get the size of the segment
   if size_seg == length(cm) 
      pos_L = find(L==L(a+1-level,1));	   
      start = L(1,1)*L(1,2) + sum(L(1:pos_L-1))*L(1,2) + 1; % always count the 1st entry twice, because that considers the size of Lo + Hi, and mult by L(1,2) becasue that size will always be the same.
      finish = (start + size_seg - 1); 
      V = [start:finish];
   else
      error('Please enter the correct level for the length of the vector.');
   end
else 
   error('Please enter the correct term for the appropriate subband. see help ICWTBAND');
end

return

%==========================================================================================
%				**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function z = c2q1d(x)

% An internal function to convert a 1D Complex vector back to a real array, 
% which is twice the height of x.
[a b] = size(x);
z = zeros(a*2,b);
skip = 1:2:(a*2);
z(skip,:) = real(x);
z(skip+1,:) = imag(x);

return