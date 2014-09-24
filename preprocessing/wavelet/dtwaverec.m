function Z = dtwaverec(C,L,biort,qshift);

% Function to perform a n-level dual-tree complex wavelet (DTCWT)
% 1D reconstruction.
%
% Z = dtwaverec(C,L,biort,qshift);
%    
%     C -> The column vector containing the Subbands
%     S -> A reference matrix containing the size of the appropriate subbands
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
%     Z -> Reconstructed signal column vector
% 
% For example:  Z = dtwaverec(C,L,'near_sym_b','qshift_b');
% performs a 3-level reconstruction from the vector C using the 13,19-tap
% near symmetric filters for level 1 and the Q-shift 14-tap filters for levels >= 2.
%
% The function ICWTBAND should be used to insert individual subbands into C
% from real or complex format vectors, before DTWAVEREC is used.
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

[a b] = size(L);

level = a; 	% No of levels = no of rows in L.

while level >= 2;  %this ensures that for level 1 we never do the following
   lh = L(a-level+1,1)*2*L(a-level+1,2); 
   if level == a,
      Lo = cwtband(C,L,level,'l','real');
   else
      Lo = Z;
   end
   Hi = cwtband(C,L,level,'h','real');
   Z = colifilt(Lo, g0b, g0a) + colifilt(Hi, g1b, g1a);
   
   if size(Z,1) ~= L(a-level+2,:)  %if Z is not the same length as the next entry in L => t1 was extended.
      Z = Z(2:size(Z,1)-1,:);     	%now we have to clip Z so that it is the same height as the next Hi
   end
   if any(size(Z) ~= L(a-level+2,:)),
      error('L matrix entries are not valid for DTWAVEREC');
   end
   
   level = level - 1;
end

if level == 1;
   if level == a,
      Lo = cwtband(C,L,1,'l','real');
   else
      Lo = Z;
   end
   Hi = cwtband(C,L,1,'h','real');
   Z = colfilter(Lo,g0o) + colfilter(Hi,g1o);
end

return


