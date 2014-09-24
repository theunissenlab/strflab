function Z = dtwaverec2(C,S,biort,qshift);

% Function to perform an n-level dual-tree complex wavelet (DTCWT)
% 1D reconstruction.
%
% Z = dtwaverec2(C,S,biort,qshift);
%    
%     C -> The column vector containing the Subbands
%     S -> A reference matrix containing the size of the appropriate subbands
%
%     biort ->  'antonini'   => Antonini 9,7 tap filters.
%               'legall'     => LeGall 5,3 tap filters.
%               'near_sym_a' => Near-Symmetric 5,7 tap filters.
%               'near_sym_b' => Near-Symmetric 13,19 tap filters.
%
%     qshift -> 'qshift_a' => Quarter Sample Shift Orthogonal Even-length (Q-Shift) 10,10 tap filters.
%               'qshift_b' => Q-Shift 14,14 tap filters.
%               'qshift_c' => Q-Shift 16,16 tap filters.
%               'qshift_d' => Q-Shift 18,18 tap filters.
%
%     Z -> Reconstructed real image matrix
%
% 
% For example:  Z = dtwaverec2(C,S,'antonini','qshift_c');
% performs a 3-level reconstruction from the vector C using the Antonini filters 
% for level 1 and the Q-shift 16-tap filters for levels >= 2.
%
% The function ICWTBAND2 should be used to insert individual subimages into C
% from real or complex format images, before DTWAVEREC2 is used.
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
      error('Please enter the correct names of the Biorthogonal or Q-Shift Filters, see help DTWAVEREC for details.');
   end
else
   error('Please enter the names of the Biorthogonal or Q-Shift Filters as shown in help DTWAVEREC.');
end

% start = 1;
finish = 0;
[a b] = size(S);

current_level = a;

while current_level >= 2; ;  %this ensures that for level -1 we never do the following
   
   sx = S(a-current_level+1,:);  %sx contains the size of the subband at the current level
   sxx = prod(sx);
   t = [1:sxx] + finish; 
   if current_level == a,
      ll = reshape(C(t),sx); 
      t = t + sxx;
   else
      ll = Z;
   end
   lh = reshape(C(t),sx); 
   t = t + sxx; 
   hl = reshape(C(t),sx);
   t = t + sxx; 
   hh = reshape(C(t),sx);
   finish = t(end);
   
   % Do even Qshift filters on columns.
   y1 = colifilt(ll,g0b,g0a) + colifilt(lh,g1b,g1a);
   y2 = colifilt(hl,g0b,g0a) + colifilt(hh,g1b,g1a);
   % Do even Qshift filters on rows.
   Z = (colifilt(y1.',g0b,g0a) + colifilt(y2.',g1b,g1a)).'; 
   
   % Check size of Z and crop as required
   [row_size col_size] = size(Z);
   if row_size ~= S(a-current_level+2,1)		%check to see if this result needs to be cropped for the rows
      Z = Z(2:row_size-1,:);
   end 
   if col_size ~= S(a-current_level+2,2)		%check to see if this result needs to be cropped for the cols
      Z = Z(:,2:col_size-1);
   end 
   if any(size(Z) ~= S(a-current_level+2,:)),
      error('S matrix entries are not valid for DTWAVEREC2');
   end
   
   current_level = current_level - 1;
end

if current_level == 1;
   
   sx = S(a-current_level+1,:);  %sx contains the size of the subband at the current level
   sxx = prod(sx);
   t = [1:sxx] + finish; 
   if current_level == a,
      ll = reshape(C(t),sx); 
      t = t + sxx;
   else
      ll = Z;
   end
   lh = reshape(C(t),sx); 
   t = t + sxx; 
   hl = reshape(C(t),sx);
   t = t + sxx; 
   hh = reshape(C(t),sx);
   finish = t(end);
   
   % Do odd top-level filters on columns.
   y1 = colfilter(ll,g0o) + colfilter(lh,g1o);
   y2 = colfilter(hl,g0o) + colfilter(hh,g1o);
   % Do odd top-level filters on rows.
   Z = (colfilter(y1.',g0o) + colfilter(y2.',g1o)).';
   
end

return

