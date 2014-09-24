function Z = dtwaverec2mw(C,S,biort,qshift);
%function Z = dtwaverec2mw(C,S,biort,qshift);
%
% Modified by MW to accept biorthogonal filters and Q-shift filters
% as input arguments. This prevent loading from disk which is slow.
%
% Function to perform an n-level dual-tree complex wavelet (DTCWT)
% 2D reconstruction.
%
% SEE ALSO: cwtInv, dtwaverec2
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
   if row_size ~= S(a-current_level+2,1)  %check to see if this result needs to be cropped for the rows
      Z = Z(2:row_size-1,:);
   end 
   if col_size ~= S(a-current_level+2,2)  %check to see if this result needs to be cropped for the cols
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

