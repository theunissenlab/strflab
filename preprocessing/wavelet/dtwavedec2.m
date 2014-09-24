function [C,S] = dtwavedec2(X,level,biort,qshift);

% Function to perform a n-level DTCWT-2D decompostion on a 2D matrix X
%
% [C,S] = dtwavedec2(X,no_levels,'Names of Biort','Name of q_shift');
%
%     X -> 2D real matrix/Image
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
%     C -> The real column vector containing the Subbands
%     S -> The "Bookkeeping" matrix containing the size of the appropriate subbands, 
%          where the last row is the orginal size of X.  The number of rows in S
%          equals the number of levels of decomposition.
% 
% Example: [C,S] = dtwavedec2(X,3,'antonini','qshift_c');
% performs a 3-level transform on the real image X using the Antonini filters 
% for level 1 and the Q-shift 16-tap filters for levels >= 2.
%
% The function CWTBAND2 should be used to extract individual subimages from C
% in real or complex format.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000


if isstr(biort) & isstr(qshift)		%Check if the inputs are strings
   biort_exist = exist([biort '.mat']);
   qshift_exist = exist([qshift '.mat']);
   if biort_exist == 2 & qshift_exist == 2;        		%Check to see if the inputs exist as .mat files
      load (biort);
      load (qshift);
   else
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEDEC2 for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEDEC2.');
end 

sx = size(X);
if any(rem(sx,2)),			%ensure that X is a even length, thus enabling it to be extended if needs be.
   error('Size of X must be a multiple of 2');
end

if level == 0, return; end

C = []; 	%initialise
S = [];

if level >= 1,
      
   % Do odd top-level filters on rows.
   Lo = colfilter(X,h0o).';
   Hi = colfilter(X,h1o).';
   
   % Do odd top-level filters on columns.
   LoLo = colfilter(Lo,h0o).';			% LoLo
   LoHi1 = colfilter(Hi,h0o).';			% LoHi => Horizontal
   HiLo1 = colfilter(Lo,h1o).';			% HiLo => Vertical
   HiHi1 = colfilter(Hi,h1o).';	      % HiHi => Diagonal  
   % C = [ LoHi1(:) ; HiLo1(:); HiHi1(:); C]; %save them in order, but do it at end for speed.
   S = [ size(LoLo) ;S];
end

if level >= 2;
   for count = 2:level;
      [row_size col_size] = size(LoLo);
      if any(rem(row_size,4)),		% Extend by 2 rows if no. of rows of LoLo are divisable by 4;
         LoLo = [LoLo(1,:); LoLo; LoLo(end,:)];
      end 
      if any(rem(col_size,4)),		% Extend by 2 cols if no. of cols of LoLo are divisable by 4;
         LoLo = [LoLo(:,1)  LoLo  LoLo(:,end)];
      end 
      
      sx = size(LoLo);
      t1 = 1:sx(1); t2 = 1:sx(2);
      sr = sx/2;
      s1 = 1:sr(1); s2 = 1:sr(2);
      
      % Do even Qshift filters on rows.
      Lo = coldfilt(LoLo,h0b,h0a).';
      Hi = coldfilt(LoLo,h1b,h1a).';
      
      % Do even Qshift filters on columns.
      LoLo = coldfilt(Lo,h0b,h0a).';	%LoLo
      LoHi = coldfilt(Hi,h0b,h0a).';	%LoHi => Horizontal
      HiLo = coldfilt(Lo,h1b,h1a).';	%HiLo => Vertical
      HiHi = coldfilt(Hi,h1b,h1a).';	%HiHi => Diagonal   
      C = [ LoHi(:) ; HiLo(:); HiHi(:); C];
      S = [ size(LoLo) ;S];
   end
end

C = [LoLo(:) ; C; LoHi1(:) ; HiLo1(:); HiHi1(:)];  % finally add in LoLo and level 1 subbands.

return

