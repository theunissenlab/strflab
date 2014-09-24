function [stim, params] = preprocSpectra(rawStim, params);
% function [stim, params] = preprocSpectra(rawStim, params);
%
% A script to preprocess movie stimuli (x-y-t matrix) into 3D Fourier spectra
%
%
% INPUT:
%   [rawStim] = A X-by-Y-by-T matrix containing stimuli (movie)
%    [params] = structure that contains parameters for preprocessing
%      .tsize = Number of frames to calculate spectra (default: 10)
% .gausslimit = Fraction of gaussian window at the edge (default: 0.01)
%    .centers = A vector ([m n]) to restrict spectra up to m-by-n low
%               spatial frequency signals.
%               (default: [0 0], no restriction)
%  .phasemode = A parameter to specify how to deal with phase information
%               0: absolute spectral amplitude (default)
%               1: spectral power
%               2: sin and cosine amplitudes
%                  (this will return 2x numbers of Fourier channels)
%               3: half-rectified signals for sin and cos amplitudes
%                  (this will return 4x numbers of Fourier channels)
%  .smoothdim = A vector to specify which dimension to be smoothed.
%               e.g., [1 1 0] will return a smoothness prior matrix A
%               as a field for PARAMS in such a way that signals
%               will be smoothed for fx and fy dimension but not
%               for ft dimension.
%               Default is NaN. (don't return a smoothness prior matrix)
%  .normalize = A flag to normalize mean and stds for each channel (default: 1)
%
% OUTPUT:
%      [stim] = Preprocessed stimuli that can be used for STRF fitting.
%               NxD matrix, N=sample size, D=dimensions.
%    [params] = structure that contains parameters for preprocessing, with additional fields:
%      .fSize = Dimensions of Fourier spectra
%      .nChan = Number of total channels
%          .A = A smoothness prior matrix. This will be only returned if .smoothdim
%               is set.
%
% EXAMPLE:
% [PARAMS] = preprocSpectra();
%   returns a structure with the default set of parameters.
%
% [stim, PARAMS] = preprocSpectra(rawStim, PARAMS)
%   preprocesses stimuli rawStim into stim using params.
%
% SEE ALSO: preprocWavelets
% ====================


%%
%% Settings for parameters
%%

if ~exist('params','var')
     params = [];
end


if ~isfield(params, 'tsize')
	 params.tsize = 10;
end


if ~isfield(params, 'gausslimit')
	 params.gausslimit = 0.05;
end


if ~isfield(params, 'centers')
	 params.centers = [0 0];
end


if ~isfield(params, 'phasemode')
    params.phasemode = 0;
end


if ~isfield(params, 'smoothdim')
    params.smoothdim = NaN;
end


if ~isfield(params, 'normalize')
	 params.normalize = 1;
end


params.class = 'preprocSpectra';

if nargin < 1  %% just called to get default parameters
    stim = params;
    return;
end



%%
%% Preprocessings
%%

ftsize = floor(params.tsize/2)+1;

moviesize = size(rawStim);

xytgauss = mymakegauss([moviesize(1:2) params.tsize], params.gausslimit);

tsize = params.tsize;
phasemode = params.phasemode;
centers = params.centers;

divdisp = ceil(moviesize(3)/20);


for ii=1:moviesize(3)
	trange = mod( (1:tsize) + ii - ceil(tsize/2)-1, moviesize(3) ) + 1;

	singled_tpatch = single(rawStim(:,:,trange));

	singled_tpatch = singled_tpatch.*xytgauss;
	
	fpatch = fftshift(fftn(singled_tpatch));
	fpatch = fpatch(:,:,end-ftsize+1:end);

	if centers
		fpatch = fpatch(ceil((end-centers(1))/2)+(1:centers(1)), ceil((end-centers(2))/2)+(1:centers(2)), :);
	end

	switch phasemode
		case 0
			fpatch = abs(fpatch);
		case 1
			fpatch = abs(fpatch).^2;
		case 2
			fpatch = cat(4, real(fpatch), imag(fpatch));
		case 3
			fpatch = cat(4, real(fpatch), imag(fpatch), -real(fpatch), -imag(fpatch));
			fpatch(find(fpatch<0)) = 0;
	end

	if ii==1
		chnum = length(fpatch(:));
		stim = zeros(moviesize(3), chnum, 'single');
	end

	stim(ii,:) = fpatch(:)';

	if mod(ii, divdisp) == 0
		fprintf('.');
	end
end

fprintf('\n');

fSize = size(fpatch);
if length(fSize) < 4
	fSize = [fSize 1];
end

params.fSize = fSize;

params.nChan = prod(fSize);

if ~any(isnan(params.smoothdim))
    params.A = getSmoothnessPrior(params.fSize, params.smoothdim);
end

if params.normalize
	if isfield(params, 'means') % Already preprocessed. Use the means and stds
		[stim] = norm_std_mean(stim, params.stds, params.means);
	else
		[stim, stds, means] = norm_std_mean(stim);
		params.means = means;
		params.stds = stds;
	end
end
