function y = cwtband(C,L,level,passband,real_or_cplx);

% Function to retrieve the subband required from the DT CWT vector C.
% y = cwtband(C,L,level,passband,real_or_cplx);
%
%     C -> The column vector containing the Cwtbands (from DTWAVEDEC.M)
%     L -> A bookkeeping matrix containing the size of the appropriate cwtbands
%     level -> The level of the required subband
%     passband ->  The type of the required subband, 'h' for Hi-pass or 'l' for Low-pass
%     **NB 'l' only works if level == size(L,1), the max level of the decomposition.
%
%     real_or_cplx ->  'real' => If y is to be a real vector
%                      'cplx' => If y is to be a complex vector. This is the DEFAULT.
%
%     y -> the output vector
%
% For example: Lo = cwtband(C,L,4,'l','real'); (only if size(L,1) = 4)
%              Hi = cwtband(C,L,2,'h');
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000

%set the default for real_or_cplx
if nargin < 5, 
   real_or_cplx = 'cplx';
end 

[a b] = size(L);

if level > (a); 	%check level is correct
   error('Error - Please Check that you are entering the correct level');
end

%get the appropiate signal vector

if passband == 'l'
   if level ~= a;
      warning(sprintf('Only the Lo cwtband of the highest level is returned, in this case level %d', a));
   end
   %reconstruct the n-x-m vector if needs be
   y = reshape(C(1:prod(L(1,:))),L(1,:));
   if real_or_cplx == 'cplx'
      t = 1:2:L(1,1);   %we know that the Lo is always going to be the first entry
      y = y(t,:)+j*y(t+1,:);
   elseif real_or_cplx == 'real'
   else
      error('Please input the correct term for REAL or COMPLEX output');
   end
   
elseif passband == 'h'
   pos_L = a+1-level;
   start = prod(L(1,:)) + sum(prod(L(1:pos_L-1,:).')) + 1; %always count the 1st entry twice, because that considers the size of Lo + Hi, and mult by L(1,2) becasue that size will always be the same.
   finish = start + prod(L(pos_L,:)) - 1;
   y = reshape(C(start:finish),L(pos_L,:));
   if real_or_cplx == 'cplx'
      t = 1:2:L(pos_L,1);
      y = y(t,:)+j*y(t+1,:);
   elseif real_or_cplx == 'real'
   else
      error('Please input the correct term for REAL or COMPLEX output');
   end
else 
   error('Please enter the correct term for the appropriate cwtband. see help CWTBAND');
   
end