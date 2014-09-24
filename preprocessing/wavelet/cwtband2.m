function Z = cwtband2(C,S,level,orientation,real_or_cplx)

% 2-D Dual-tree Complex Wavelet Transform:
% Function to retrieve the subimage required from the 2-D DT CWT vector C.
% 
% output = cwtband2(C,S,level,orientation,real_or_cplx)
%
%     C -> The column vector containing the Subbands
%     S -> The "Bookkeeping" matrix
%     level -> The level which you want to get a subimage from
%
%     orientation -> 'l' => LoLo subimage (Only at the highest transformed level and can return a complex matrix)
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
% Example:     ll = cwtband2(C,S,4,'l','real');
%              lh = cwtband2(C,S,2,'h','cplx');
%              hl = cwtband2(C,S,1,'v','cplx');
%              hh = cwtband2(C,S,3,'d','real');
%
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

if orientation == 'l';   
   if level ~= a  %General warning
      warning(sprintf('Only the LoLo subband of the highest level is returned, in this case level %d', a));
   end
   Z = reshape(C(1:prod(S(1,:))),S(1,:));
   if real_or_cplx == 'real'			
   elseif real_or_cplx == 'cplx'	
      Z = q2c(Z);
   else
      error('Unsupported reference when calling CWTBAND2, relating to REAL or CPLX')
   end
   return
   
else orientation ~= 'l'; %for all subimages except LoLo
   ysize = S(a+1-level,:);
   row_size = ysize(1);
   col_size = ysize(2);
   size_seg = row_size*col_size;   							%get the size of the segment
   row_pos = find(S(:,1)==[row_size]);	     %the position of the of the row_size in S
   col_pos = find(S(:,2)==[col_size]);	 
   start_sec = S(1,1)*S(1,2) + sum(S(1:row_pos-1,1).*S(1:col_pos-1,2))*3; %Now start is the start of chosen section
end

if orientation == 'h' %| orientation == 'h2';
   start = start_sec + 1; %Because for Horziontal its is the first one in the list
   finish = (start + size_seg - 1); 
   Z = reshape(C(start:finish),ysize);
   if real_or_cplx == 'cplx'
      Z = q2c(Z);
   elseif real_or_cplx == 'real'		
   else
      error('Unsupported reference when calling CWTBAND2, relating to REAL or CPLX')
   end
   
elseif orientation == 'v' %| orientation == 'v2';
   start = start_sec + size_seg + 1;
   finish = (start + size_seg - 1); 
   Z = reshape(C(start:finish),ysize);
   if real_or_cplx == 'cplx'
      Z = q2c(Z);
   elseif real_or_cplx == 'real'				
   else
      error('Unsupported reference when calling CWTBAND2, relating to REAL or CPLX')
   end
   
elseif orientation == 'd' %| orientation == 'd2';
   start = start_sec + size_seg*2 + 1;
   finish = (start + size_seg - 1); 
   Z = reshape(C(start:finish),ysize);
   if real_or_cplx == 'cplx'
      Z = q2c(Z);
   elseif real_or_cplx == 'real'					
   else
      error('Unsupported reference when calling CWTBAND2, relating to REAL or CPLX')
   end    
else
   error('Unsupported reference when calling CWTBAND2, relating to Orientation')
end
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
z = [a-d; a+d]*sqrt(0.5) + [b+c; b-c]*sqrt(-0.5);

return
