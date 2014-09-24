function A = getSmoothnessPrior(hsize, sdimension);
%function A = getSmoothnessPrior(hsize, sdimension)
%
% A function to make a smoothnessprior matrix for N-D matrix of size: hsize
% The matrix can be either 1D ~ 3D or 4D
%
% INPUT:
%      [hsize] = vector of sizes for each dimension of matrix
% [sdimension] = determines which dimensions to be smoothed
%
% OUTPUT:
%          [A] = sparse smoothness prior matrix
%
% EXAMPLE:
%   A = getSmoothnessPrior([20 20 5], [1 1 0]);
%


if length(hsize) <= 4
  hsize = [hsize ones(1,4-length(hsize))];
end


% The first assumption is to smooth all 1st to 3rd dimensions

smdim = 3;  % 3x3x3 smoother matrix
s = zeros(smdim, smdim, smdim);
s(:,:,1) = [0 0 0; 0 -1 0; 0 0 0];
s(:,:,2) = [0 -1 0; -1 6 -1; 0 -1 0];
s(:,:,3) = [0 0 0; 0 -1 0; 0 0 0];
scent = ceil(smdim/2);
ssize = size(s);



if exist('sdimension') % smooth only specified dimenstions
	for sdir=1:3
		if sdimension(sdir)
			ssubs{sdir} = 1:smdim;  % smooth for this dimension
		else
			ssubs{sdir} = scent;  % don't smooth for this dimension
		end
	end
	s2 = zeros(smdim, smdim, smdim);
	s2(ssubs{1}, ssubs{2}, ssubs{3}) = s(ssubs{1}, ssubs{2}, ssubs{3});

	s2(scent, scent, scent) = 0;
	s = s2;
	s(scent, scent, scent) = -sum(s(:));
end

% Now the smoother matrix was modified to smooth only specified dimension(s)

nonzeroind = find(s);
nonzeroval = s(nonzeroind);
numnonzeros = length(nonzeroind);

for jj=1:numnonzeros
  [sind(jj,1) sind(jj,2) sind(jj,3)] = ind2sub(ssize, nonzeroind(jj));
end
sind = sind - scent;

dimension3 = prod(hsize(1:3));
dimension4 = prod(hsize);

hsize3 = hsize(1:3);

% make a smoothness prior matrix
%A = zeros(dimension4, dimension4, 'single');
A = sparse([],[],[],dimension4, dimension4, dimension4*length(find(s))*hsize(4));
for ii=1:dimension3
  [d(1) d(2) d(3)] = ind2sub(hsize(1:3), ii);

  for jj=1:numnonzeros
    thisind = d+sind(jj,:);
    if checkwithin(hsize3, thisind)
      ind = sub2ind(hsize, thisind(1), thisind(2), thisind(3));
      A(ii,ind) = nonzeroval(jj);
    end
  end
end


%% expand to the 4th dimension
if hsize(4)>=2
  for ii=2:hsize(4)
    thisdim = (1:dimension3)+dimension3*(ii-1);
    A(thisdim, thisdim) = A(1:dimension3, 1:dimension3);
  end
end


%% sum for each dimension should be 0.
asum = sum(A,2);
diags = find(diag(ones(1,size(A,1))));
A(diags) = A(diags) - asum;


return;


%%%-----------------------------
function out = checkwithin(hsize, indarray)

if any(find(indarray<=0))
  out = 0;
  return;
end

indarray2 = indarray - hsize;
if any(find(indarray2>0))
  out = 0;
  return;
end

out = 1;
