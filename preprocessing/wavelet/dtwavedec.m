function [C,L] = dtwavedec(X,level,biort,qshift);

% Function to perform a n-level dual-tree complex wavelet (DTCWT)
% decomposition on a 1D column vector X.  The filter names are
% given by biort and qshift.  The length of X must be even
% but does NOT need to be a power of 2.
%
% [C,L] = dtwavedec(X,level,biort,qshift);
%
%     X -> Signal Column vector
%
%     level -> No. of levels of wavelet decomposition
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%
%     qshift -> 'qshift_a' => Quarter Sample Shift Orthogonal Even-length 
%                             (Q-Shift) 10,10 tap filters (only 6,6 non-zero taps).
%               'qshift_b' => Q-Shift 14,14 tap filters.
%               'qshift_c' => Q-Shift 16,16 tap filters.
%               'qshift_d' => Q-Shift 18,18 tap filters.
%
%     C -> The column vector containing the Subbands
%     L -> The "Bookkeeping" matrix containing the size of the appropriate subbands, 
%          where the last row is the orginal size of X.  The number of rows in L
%          equals the number of levels of decomposition.
%
% For example:  [C,L] = dtwavedec(X,3,'near_sym_b','qshift_b');
% performs a 3-level transform on the real column vector X using the 13,19-tap
% near symmetric filters for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% The function CWTBAND should be used to extract individual subbands from C
% in real or complex format.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000


if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2;  %Check to see if the filters exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEDEC for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEDEC.');
end

C = [];			% initialise
L = size(X);

if any(rem(L(1),2)),	 % ensure that X is an even length, thus enabling it to be extended if needs be.
   error('Size of X must be a multiple of 2');
end

if level == 0, return; end

if level >= 1;   %always have to do level 1  
   Hi = colfilter(X, h1o);   
   Lo = colfilter(X, h0o);
   C = [Hi(:) ; C];
end

if level >= 2;
 
   for count = 2:level;  
      if rem(size(Lo,1),4),		%check to see if Lo is divisable by 4, if not extend.
         Lo = [Lo(1,:); Lo; Lo(end,:)];
      end     
      
      sx = size(Lo);    %sx = size of Lo
      L = [ sx(1)/2 sx(2); L];  
      Hi = coldfilt(Lo,h1b,h1a);
      Lo = coldfilt(Lo,h0b,h0a); 
      C = [Hi(:) ; C];      
   end   
end
C = [Lo(:); C];

return
