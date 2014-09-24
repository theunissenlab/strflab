function out = ndimages(im, Params)

%function out = ndimages(im, Params)
%
% Allows for the display of n-dimensional (up to 6-d) images by tiling
% along higher dimensions
%
% INPUT:
%      [im] = a matrix to be displayed (up to 6 dimensions)
%  [Params] = structure that contains parameters
%     .clim = The minimum and maximum value for color plots. e.g., [0 1]
%    .title = A string to show on top of the images
% .X and .Y = Strings to show as labels for large axes
% .x and .y = Strings to show as labels for small axes
%  .fftgrid = A flag to add zero frequency axes (assuming FFT-ed data)
%    .imdim = Force the matrix dimension as this value
%
% OUTPUT:
%     [out] = outputs first tiled image
%
% EXAMPLE:
%   ndimages(IM)
%   ndimages(IM, [CLOW CHIGH])
%   ndimages(IM, PARAMS)
%

imsize = size(im);

fftgrid = 0;
axisimage = 0;
figshow = 1;
axisxy = 0;

if ~exist('Params', 'var')
	cmin = min(im(:));
	cmax = max(im(:));
elseif ~isstruct(Params)
	cmin = Params(1);
	cmax = Params(2);
else
	%Params is actually a structure of parameters
	if isfield(Params, 'clim')
		cmin = Params.clim(1);
		cmax = Params.clim(2);
	else
		cmin = min(im(:));
		cmax = max(im(:));
	end
	if isfield(Params, 'title')
		ftitle = Params.title;
	end
	if isfield(Params, 'X')
		label_lx = Params.X;
	end
	if isfield(Params, 'Y')
		label_ly = Params.Y;
	end
	if isfield(Params, 'x')
		label_sx = Params.x;
	end
	if isfield(Params, 'y')
		label_sy = Params.y;
	end
	if isfield(Params, 'axisimage')
		axisimage = Params.axisimage;
	end
	if isfield(Params, 'fftgrid')
		fftgrid = Params.fftgrid;
	end
	if isfield(Params, 'imdim')
		if length(imsize) < Params.imdim
			imsize = [imsize ones(1,Params.imdim-length(imsize))];
		end
	end
	if isfield(Params, 'show')
		figshow = Params.show;
	end
	if isfield(Params, 'axisxy')
        axisxy = Params.axisxy;
	end
end

if (cmin==0 & cmax==0) | cmin==cmax
	cmin=0;
	cmax=1;
end

if figshow, clf, end

if length(imsize) == 2
	if figshow
		imagesc(im, [cmin cmax]);
	end
	out = im;

elseif length(imsize) == 3 | length(imsize) == 4
	if length(imsize) == 3
		imnum = imsize(3);
		sx = ceil(imnum^0.5);
		sy = ceil(imnum/sx);
	else
		imnum = prod(imsize(3:4));
		sx = imsize(3);
		sy = imsize(4);
		im = reshape(im, [imsize(1) imsize(2) imsize(3)*imsize(4)]);
	end
	for ii=1:imnum
		subplot(sy, sx, ii);
		imagesc(im(:,:,ii), [cmin cmax]);
		settickoff();
		if fftgrid
			sz = size(im);
			hold on;
			h = plot([0 sz(1)+1], [ceil(sz(2)/2) ceil(sz(2)/2)], 'r--');
			h2 = plot([ceil(sz(1)/2) ceil(sz(1)/2)], [0 sz(2)], 'r--');
			if length(fftgrid) == 3
				set(h, 'color', fftgrid);
				set(h2, 'color', fftgrid);
			end
		end
		if ii == ceil(sx/2) & exist('ftitle', 'var')
			title(ftitle);
			settickoff();
		end
		if ii == sx*(sy-1)+1
			if exist('label_sx', 'var'), xlabel(label_sx), end
			if exist('label_sy', 'var'), ylabel(label_sy), end
			settickoff();
		end
		if ii == round(sy/2)*sx+1 & exist('label_ly', 'var')
			h = ylabel(label_ly);
			set(h, 'fontweight','bold');
			settickoff();
		end
		if ii == sy*sx-round(sx/2)+1 & exist('label_lx', 'var')
			h = xlabel(label_lx);
			set(h, 'fontweight','bold');
			settickoff();
		end

		if axisimage, axis image, end
		if axisxy, axis xy, end
	end
else % 5 or 6 dimensions
	sx = imsize(3);
	sy = imsize(4);
	bigims = ones(sy*imsize(1)+sy, sx*imsize(2)+sx)*cmax;
	if length(imsize)==5
		imnum = imsize(5);
		lx = ceil(imnum^0.5);
		ly = ceil(imnum/lx);
	else
		imnum = imsize(5)*imsize(6);
		lx = imsize(5);
		ly = imsize(6);
		im = reshape(im, [imsize(1) imsize(2) imsize(3) imsize(4) imsize(5)*imsize(6)]);
	end
	for lcount=1:imnum
		jx = mod(lcount-1, lx)+1;
		jy = floor((lcount-1)/lx)+1;
		for ty=1:imsize(4)
			for tx=1:imsize(3)
				xrange = (1:imsize(2))+(tx-1)*(1+imsize(2));
				yrange = (1:imsize(1))+(ty-1)*(1+imsize(1));
				bigims(yrange,xrange) = im(:,:,tx, ty, lcount);
			end
		end
		if figshow
			subplot(ly, lx, lcount);
			imagesc(bigims, [cmin cmax]);
			if axisimage, axis image, end
			if axisxy, axis xy, end
			settickoff();
			if jx==1 & jy==prod(imsize(5))
				if exist('label_sx')
					xlabel(label_sx);
				end
				if exist('label_sy')
					ylabel(label_sy);
				end
			elseif jx==1 & jy==ceil(ly/2) & exist('label_ly', 'var')
				h = ylabel(label_ly);
				set(h, 'FontSize', 20);
			elseif jx==ceil(lx/2) & jy==prod(ly) & exist('label_lx', 'var')
				h = xlabel(label_lx);
				set(h, 'FontSize', 20);
			elseif jx==ceil(lx/2) & jy==1 & exist('ftitle', 'var')
				title(ftitle);
			end
		end
	end
end

if exist('bigims', 'var')
	out = bigims;
end

if ~exist('out','var')
	out = [];
end

return

%%%%%%%--------------------------------
function settickoff()

set(gca, 'xtick', []);
set(gca, 'ytick', []);
axis on;

return
