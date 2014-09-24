function Z = cwtband6(C,S,level)

% 2-D Dual-tree Complex Wavelet Transform:
% Function to retrieve the set of 6 complex bandpass subimages 
% at level from the 2-D DT CWT vector C.
% 
% output = cwtband6(C,S,level)
%
%     C -> The column vector containing the Subbands
%     S -> The "Bookkeeping" matrix
%     level -> The level the submages are to come from.
%
% The 6 subimages are returned in Z(:,:,1:6) in the order:
%    h1, h2, v1, v2, d1, d2 
% (oriented at +15, -15, +75, -75, +45, -45 degrees).
%
% Example: 
%
%    bp = cwtband6(C,S,2);
%    figure;cimage5([bp(:,:,1) bp(:,:,3) bp(:,:,5);bp(:,:,2) bp(:,:,4) bp(:,:,6)])
%
% returns the 6 complex subimages at level 2 in bp(:,:,1:6) and plots them.
%
% Cian Shaffrey, Nick Kingsbury
% Cambridge University, Sept 2000

[a,b] = size(S);

if level > (a); 	%check level is correct
   error('Error - Please Check that you are entering the correct level');
end


ysize = S(a+1-level,:);
row_size = ysize(1);
col_size = ysize(2);
size_seg = row_size*col_size;   							%get the size of the segment
row_pos = find(S(:,1)==[row_size]);	     %the position of the of the row_size in S
col_pos = find(S(:,2)==[col_size]);	 
start_sec = S(1,1)*S(1,2) + sum(S(1:row_pos-1,1).*S(1:col_pos-1,2))*3; %Now start is the start of chosen section

Z = zeros(row_size/2,col_size/2,6);

% Horizontal band pair:
start = start_sec + 1; %Because for Horziontal its is the first one in the list
finish = (start + size_seg - 1); 
Z(:,:,1:2) = q2c(reshape(C(start:finish),ysize));

% Vertical band pair:
start = start_sec + size_seg + 1;
finish = (start + size_seg - 1); 
Z(:,:,3:4) = q2c(reshape(C(start:finish),ysize));

% Diagonal band pair:
start = start_sec + size_seg*2 + 1;
finish = (start + size_seg - 1); 
Z(:,:,5:6) = q2c(reshape(C(start:finish),ysize));

return

%==========================================================================================
%						**********  	INTERNAL FUNCTION    **********
%==========================================================================================

function z = q2c(y)

% function z = q2c(y)
% Convert from quads in y to complex numbers in z.

sy = size(y);
t1 = 1:2:sy(1); t2 = 1:2:sy(2);

% Arrange pixels from the corners of the quads into
% 2 subimages of alternate real and imag pixels.
%  a----b
%  |    |
%  |    |
%  c----d
a = y(t1,t2);
b = y(t1,t2+1);
c = y(t1+1,t2);
d = y(t1+1,t2+1);

% Form the real and imag parts of the 2 subbands.
z(:,:,1) = [a-d]*sqrt(0.5) + [b+c]*sqrt(-0.5);
z(:,:,2) = [a+d]*sqrt(0.5) + [b-c]*sqrt(-0.5);

return
