function [C,S] = dtwavedec2mw(X,level,biort,qshift);
%function [C,S] = dtwavedec2mw(X,level,biort,qshift);
%
% Modified by MW to accept biorthogonal filters and Q-shift filters
% as input arguments. This prevent loading from disk which is slow.
%
% Function to perform a n-level DTCWT-2D decompostion on a 2D matrix X
%
% SEE ALSO: cwtDT, dtwavedec2
%
% By Michael Wu  --  waftingpetal@yahoo.com (Mar 2007)
%
% ====================


% Init Antonini & qShiftC filters
%--------------------
g0o=biort.g0o;
g1o=biort.g1o;
h0o=biort.h0o;
h1o=biort.h1o;

g0a=qshift.g0a;
g0b=qshift.g0b;
g1a=qshift.g1a;
g1b=qshift.g1b;
h0a=qshift.h0a;
h0b=qshift.h0b;
h1a=qshift.h1a;
h1b=qshift.h1b;



sx = size(X);
if any(rem(sx,2)),  %ensure that X is a even length, thus enabling it to be extended if needs be.
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

