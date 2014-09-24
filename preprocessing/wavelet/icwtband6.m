function [cm,V] = icwtband6(ym,S,level);

% 2-D Dual-tree Complex Wavelet Transform:
% Function to determine where to insert a modified set of 6 subimages ym 
% into C and to convert from complex to real format.
% It is necessary to use 'C(V) = cm;' to then do the insertion.
% (For large arrays C, this is much more efficient than copying
% the whole of C into and out of the function.)
%
% [cm,V] = icwtband6(ym,S,level);
%
%     ym -> The modified subimage matrix
%     S -> The "Bookkeeping" matrix
%     level -> The level the submages are to go to.
%
%     cm -> Is the returned subimage (in vector form) ready for subistution into C according to V
%     V -> Is the vector of points where cm should be inserted in to C.
%
% For example:  [cm,V] = icwtband6(ym,S,2);
%               C(V) = cm;  
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, September 2000

[a,b] = size(S);

if level > (a); 	%check level is correct
   error('Error - Please Check that you are entering the correct level');
end

% Check the size and type of ym:
[c d e] = size(ym);
if (e ~= 6) | isreal(ym),
   error('ym must be a 3-D array of 6 complex subimages.');
end

row_size = S(a+1-level,1);
col_size = S(a+1-level,2);
if (row_size/2 ~= c) | (col_size/2 ~= d),  % 1/2 is due to the conversion of reals to cplx
   error('ym is not the correct size for the chosen level');
end

% Convert complex subimage pairs to real subimages:
yh = c2q(ym(:,:,1:2)); % Horizontal pair
yv = c2q(ym(:,:,3:4)); % Vertical pair
yd = c2q(ym(:,:,5:6)); % Diagonal pair
cm = [yh(:); yv(:); yd(:)]; % Assemble these as the output vector.

% Now to get the start and finish points
row_pos = find(S(:,1)==[row_size]);	   % the position of the of the row_size in S
col_pos = find(S(:,2)==[col_size]);	 
start = S(1,1)*S(1,2) + sum(S(1:row_pos-1,1).*S(1:col_pos-1,2))*3 + 1;
size_seg = length(cm);   					% get the size of the segment
finish = start + size_seg - 1; 
V = [start:finish];

return;

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function x = c2q(w)

% function z = c2q(w)
% Convert from complex w(:,:,1:2) to real quad-numbers in z.

% Arrange pixels from the real and imag parts of the 2 subbands
% into 4 separate subimages .
%  A----B     Re   Im
%  |    |
%  |    |
%  C----D     Re   Im

A = real(w(:,:,1));
B = imag(w(:,:,1));
C = real(w(:,:,2));
D = imag(w(:,:,2));

t1 = 1:2:(2*size(w,1));
t2 = 1:2:(2*size(w,2));

% Recover each of the 4 corners of the quads.
sc = sqrt(0.5);

x(t1,t2)     = (A+C)*sc; % a 
x(t1,t2+1)   = (B+D)*sc; % b 
x(t1+1,t2)   = (B-D)*sc; % c
x(t1+1,t2+1) = (C-A)*sc; % d

return