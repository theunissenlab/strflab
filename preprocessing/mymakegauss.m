function mygauss = mymakegauss(xytsize, sratio, showfig)
%function mygauss = mymakegauss(xytsize, sratio, showfig)
%
% Creates a 3-D gaussian
%
% INPUT:
% [xytsize] = size of a matrix to make gaussian
%  [sratio] = the maximum value at the edge of the matrix
%             determines width of gaussian
% [showfig] = whether or not to show figure of the gaussian
%
% OUTPUT:
% [mygauss] = 3-d gaussian
%

dx = -1:(2/(xytsize(1)-1)):1;
dy = -1:(2/(xytsize(2)-1)):1;

if length(xytsize) == 2
	xytsize = [xytsize 1];
end

if xytsize(3) == 1
	dt = 0;
else
	dt = -1:(2/(xytsize(3)-1)):1;
end

sigma = sqrt(-0.5/log(sratio));

[ix iy it] = ndgrid(dx, dy, dt);

mygauss = exp( - ((ix.^2)+(iy.^2)+(it.^2))/(2*sigma^2) );

if ~exist('showfig','var')
	return
end

if showfig
	clf;
	sp = ceil(sqrt(xytsize(3)));
	for ii=1:xytsize(3)
		subplot(sp, sp, ii);
		imagesc(mygauss(:,:,ii), [0 1]);
	end
end
