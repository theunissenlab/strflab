function [cm,V] = icwtband2(ym,S,level,orientation,real_or_cplx);

% 2-D Dual-tree Complex Wavelet Transform:
% Function to determine where to insert a modified subimage ym in C
% and to convert from complex to real format if ym is complex.
% It is necessary to use 'C(V) = cm;' to then do the insertion.
% (For large arrays C, this is much more efficient than copying
% the whole of C into and out of the function.)
%
% [cm,V] = icwtband2(ym,S,level,orientation,real_or_cplx);
%
%     ym -> The modified subimage matrix
%     S -> The "Bookkeeping" matrix
%     level -> The level the submages are to go to.
%     orientation -> 'l' => LoLo subimage (Only at the highest transformed level is available)
%                    'h' => Horziontal Edge Image @ +15 on top
%                        => Horziontal Edge Image @ -15 on the bottom
%                    'v' => Vertical Edge Image @ +75 on top
%                        => Vertical Edge Image @ -75 on the bottom
%                    'd' => Diagonal Edge Image @ +45 on top
%                        => Diagonal Edge Image @ -45 on the bottom
%
%     real_or_cplx-> 'real' => Return the quad-number subimage as a purely REAL matrix
%                              (In this case the images at +/- 15 etc are combined 
%                              into a single subimage of quad numbers.)
%                    'cplx' => Return the subimage pair as a COMPLEX matrix
%
%     cm -> Is the returned subimage (in vector form) ready for subistution into C according to V
%     V -> Is the vector of points where cm should be inserted in to C.
%
% For example:  [cm,V] = icwtband2(ym,S,2,'l','real');
% or:           [cm,V] = icwtband2(ym,S,1,'d','cplx');
%               C(V) = cm;  
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, August 2000

%set the default for real_or_cplx
if nargin < 5, 
   real_or_cplx = 'cplx';
end 

[a,b] = size(S);

if level > (a); 	%check level is correct
   error('Error - Please Check that you are entering the correct level');
end

%check the size
[c d] = size(ym);
row_size = S(a+1-level,1);
col_size = S(a+1-level,2);

if real_or_cplx == 'cplx';   
   if row_size == c & col_size/2 == d;  %this 1/2 is due to the conversion of reals to cplx
      ym = c2q(ym);
   else
      error('Please enter the correct level ');
   end
   
elseif real_or_cplx == 'real';
   if row_size == c & col_size == d;
   else
      error('The level and the size of ym do not agree.');
   end   
else
   error('Please enter the correct term for the appropriate cwtband. see help ICWTBAND');
end

% now ym is now a purely real array
cm = ym(:);

%now to get the start and finish points
size_seg = row_size*col_size;   							%get the size of the segment
row_pos = find(S(:,1)==[row_size]);	     %the position of the of the row_size in S
col_pos = find(S(:,2)==[col_size]);	 
start_sec = S(1,1)*S(1,2) + sum(S(1:row_pos-1,1).*S(1:col_pos-1,2))*3;

if orientation == 'l';
   if level ~= a;
      error(sprintf('Only the Lo cwtband of the highest level is to be inserted, in this case level %d', a));
   end
   V = [1:row_size*col_size];
elseif orientation == 'h';
   start = start_sec + 1;
   finish = (start + size_seg - 1); 
   V = [start:finish];
elseif orientation == 'v';
   start = start_sec + size_seg + 1;
   finish = (start + size_seg - 1); 
   V = [start:finish];
elseif orientation == 'd';
   start = start_sec + size_seg*2 + 1;
   finish = (start + size_seg - 1); 
   V = [start:finish];
else
   error('Unsupported reference when calling CWTBAND2, relating to Orientation')
end



%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function x = c2q(w)

% function z = c2q(w)
% Convert from complex w to real quad-numbers in z.

% Arrange pixels from the real and imag parts of the 2 subbands
% into 4 separate subimages .
%  A----B     Re   Im
%  |    |
%  |    |
%  C----D     Re   Im

t1 = 1:(0.5*size(w,1));
A = real(w(t1,:));
B = imag(w(t1,:));
t1 = t1 + 0.5*size(w,1);
C = real(w(t1,:));
D = imag(w(t1,:));

t1 = 1:2:size(w,1);
t2 = 1:2:(2*size(w,2));

% Recover each of the 4 corners of the quads.
sc = sqrt(0.5);

x(t1,t2)     = (A+C)*sc; % a 
x(t1,t2+1)   = (B+D)*sc; % b 
x(t1+1,t2)   = (B-D)*sc; % c
x(t1+1,t2+1) = (C-A)*sc; % d

return