function [PS, params] = preprocWavelets(S, params);
% function [PS, params] = preprocWavelets(S, params);
%
% A script for preprocessing of stimuli using a Gabor wavelet bais set
%
%  PARAMS = preprocWavelets;
%    returns the default set of parameters.
%
%  [PS, PARAMS] = preprocWavelets(S, PARAMS)
%    returns preprocessed stimuli (or wavelets) and parameters
%
%  PARAMS is a structure that contains parameters for preprocessing
%      .dirdivisions     Number of directions for wavelets (default: 8)
%      .fdivisions     Number of frequencies (default: 5)
%      .veldivisions     Number of velocities (default: 5)
%      .tsize   Number of frames to calculate wavelets (default: 10)
%      .sfmax   The maximum spatial frequency/stimulus size at zero velocity (default: 9)
%      .sfmin   The minimum spatial frequency/stimulus size  at zero velocity (default: 2)
%      .f_step_log   A flag to specify linear or log step of frequency (default: 0)
%      .tfmax   The maximum temporal frequency/stimulus size (default: 3.5)
%      .sf_gaussratio   The ratio between the Gaussian window and spatial frequency
%                       (default: 0.5)
%      .tf_gaussratio   The ratio between the Gaussian window and temporal frequency
%                       (default: 0.4)
%      .std_step   Spatial separation of each wavelet in terms of sigma of
%                  the Gaussian window (default: 2.5)
%      .phasemode   A parameter to specify how to deal with phase information
%               0: spectral amplitude
%               1: linear sin and cos phase ampliture (2x number of wavelets)
%               2: half-rectified sin and cos phase amplitude (4x number of wavelets)
%               3: 0+1 (3x number of wavelets)
%               4: 0+2 (5x number of wavelets)
%                (default: 0)
%      .phasemode_sfmax   The maximum spatial frequency to use phase information
%                         For higher frequency wavelets over this value, only
%                         spectral amplitude information are used.  (default: Inf)
%      .show_or_preprocess   If this is set to 0, the function returns wavelets
%                            of size(S) * number of channels, instead of preprocessed
%                            stimuli. This may be used for visualization purpose.
%                            If .valid_w_index is also set, this returns only a subset of
%                            wavelets specified by .valid_w_index. (default: 1)
%      .senv_max   The maximum spatial envelope (default: 0.3)
%      .tenv_max   The maximum temporal envelope (default: 0.3)
%      .local_dc   A flag to add localized dc (i.e., 0 spatial freq.) channels (default: 0)
%      .valid_w_index   This is used to specify a subset of wavelets to obtain.
%                       (See .show_or_preprocess)
%      .normalize   A flag to normalize mean and stds for each channel (default: 1)
%      .verbose    a flag for verbose mode (default: 1)
%
%  (PARAMS fileds to be modified)
%
%      .N            Number of preprocessed channels
%      .gaborparams  A set of parameters for each Gabor wavelet.
%                    This is a p-by-N matrix where p is number of parameters (8)
%                    and N is the number of wavelet channels
%                    Each field in gaborparams represents:
%                    [pos_x pos_y direction s_freq t_freq s_size t_size phasevalue]
%                    phasevalue can be 0 to 6, where
%                         0: spectra
%                         1: linear sin transform
%                         2: linear cos transform
%                         3: half-rectified sin transform (positive values)
%                         4: half-rectified sin transform (negative values)
%                         5: half-rectified cos transform (positive values)
%                         6: half-rectified cos transform (negative values)
%

%%% Set up parameters
if ~exist('params','var')
     params = [];
end
params = setDefaultParameters(params);

if nargin<1 % just called to set the default set of parameters
	PS = params;
	return;
end


start_t = cputime;

stimxytsize = size(S);
if length(stimxytsize) == 2
	stimxytsize = [stimxytsize 1]; % make sure 3 dimensions
end

patchxytsize = [stimxytsize(1:2) params.tsize];
xypixels = prod(patchxytsize(1:2));
verbose = params.verbose;

S = single(S);
S = reshape(S, [prod(stimxytsize(1:2)) stimxytsize(3)]);


%%% Make a list of gabor parameters
if verbose, fprintf('Making a list of gabor parameters... '); end

[gparams] = getGaborParameters(params);
waveletchannelnum = size(gparams,2);

if verbose, fprintf('channel num: %d\n', waveletchannelnum); end

if verbose & any(params.valid_w_index)
	fprintf('Valid channel num: %d\n', length(params.valid_w_index));
end


%%% Set up a matrix to fill-in
if params.show_or_preprocess
	if verbose, disp('Preprocessing...'); end
	PS = zeros(stimxytsize(3), waveletchannelnum, 'single');
else
	if verbose, disp('Making wavelets...'); end
	if ~any(params.valid_w_index)
		gnum = length(waveletchannelnum);
	else
		gnum = length(params.valid_w_index);
	end
	gaborbank = zeros([patchxytsize gnum], 'single');
end


%%% Preprocessing...

% ignore wavelet pixels for speed-up where:
masklimit = 0.001;   %% pixel value < masklimit AND
maskenv_below = 0.1; % spatial envelope < maskenv_below x stimulus size

lastgparam = zeros(8,1);
wcount = 0;
for ii=1:waveletchannelnum

	if any(params.valid_w_index) & ~any(ii==params.valid_w_index), continue, end

	thisgparam = gparams(:,ii);
	thesame = 1;
	if any(thisgparam(1:7) ~=lastgparam(1:7))
		thesame = 0;
	end
	if ~thesame
		[gabor0 gabor90] = make3dgabor(patchxytsize, [thisgparam(1:end-1); 0]);
		lastgparam = thisgparam;
	end
	phaseparam = thisgparam(8);
	if params.show_or_preprocess
		if ~thesame
			gabor0 = reshape(gabor0,[xypixels params.tsize]);
			gabor90 = reshape(gabor90,[xypixels params.tsize]);
			senv = thisgparam(6);
			if senv<maskenv_below
				g0 = find(max(abs(gabor0),[],2)>masklimit);
				chout0 = dotdelay(gabor0(g0,:), S(g0,:));
				g90 = find(max(abs(gabor90),[],2)>masklimit);
				chout90 = dotdelay(gabor90(g90,:), S(g90,:));
			else
				chout0 = dotdelay(gabor0, S);
				chout90 = dotdelay(gabor90, S);
			end
		end
		switch phaseparam
		case 0
			chout = sqrt(chout0.^2 + chout90.^2);
			PS(:,ii) = chout;
		case 1
			chout = chout0;
			PS(:,ii) = chout;
		case 2
			chout = chout90;
			PS(:,ii) = chout;
		case 3
			chout = chout0;
			chout(find(chout<0)) = 0;
			PS(:,ii) = chout;
		case 4
			chout = chout0;
			chout(find(chout>0)) = 0;
			PS(:,ii) = -chout;
		case 5
			chout = chout90;
			chout(find(chout<0)) = 0;
			PS(:,ii) = chout;
		case 6
			chout = chout90;
			chout(find(chout>0)) = 0;
			PS(:,ii) = -chout;
		end
	else
		wcount = wcount + 1;
		switch phaseparam
		case {0,1,3,4}
			gaborbank(:,:,:,wcount) = gabor0;
		case {2,5,6}
			gaborbank(:,:,:,wcount) = gabor90;
		end
	end

	if verbose
		if mod(ii, 50) == 0, fprintf('.'), end
		if mod(ii, 1000) == 0, fprintf(' %d channels done.\n', ii), end
	end
end
if verbose & params.show_or_preprocess, fprintf(' %d channels done.\n', ii); end


if verbose
	disp(sprintf('Wavelet preprocessing done in %.1f min.', (cputime-start_t)/60));
	if params.show_or_preprocess
		disp(sprintf('%d channels, %d samples', size(PS,2), size(PS,1)));
	else
		disp(sprintf('%d channels', size(gaborbank,4)));
	end
end


if params.show_or_preprocess
	if params.normalize
		if isfield(params, 'means') % Already preprocessed. Use the means and stds
			[PS] = norm_std_mean(PS, params.stds, params.means);
		else
			[PS, stds, means] = norm_std_mean(PS);
			params.means = means;
			params.stds = stds;
		end
	end
else % return gabors, not pre-processed data
	PS = gaborbank;
end

params.gaborparams = gparams;
params.N = size(gparams,2);

return;



%---------------------------------------------------------------------
%  Default parameter settings
%---------------------------------------------------------------------

function params = setDefaultParameters(params)

if ~isfield(params, 'dirdivisions')
	 params.dirdivisions = 8;
end

if ~isfield(params, 'fdivisions')
	 params.fdivisions = 5;
end

if ~isfield(params, 'veldivisions')
	 params.veldivisions = 5;
end

if ~isfield(params, 'tsize')
	 params.tsize = 9;
end

if ~isfield(params, 'sfmax')
	 params.sfmax = 9.0;
end

if ~isfield(params, 'sfmin')
	 params.sfmin = 2.0;
end

if ~isfield(params, 'f_step_log')
	 params.f_step_log = 0;
end

if ~isfield(params, 'tfmax')
	 params.tfmax = 3.0;
end

if ~isfield(params, 'sf_gaussratio')
	 params.sf_gaussratio = 0.5;
end

if ~isfield(params, 'tf_gaussratio')
	 params.tf_gaussratio = 0.4;
end

if ~isfield(params, 'std_step')
	 params.std_step = 2.5;
end

if ~isfield(params, 'phasemode')
	params.phasemode = 0;
	params.phasemode_sfmax = Inf;
end

if ~isfield(params, 'senv_max')
	 params.senv_max = 0.3;
end

if ~isfield(params, 'tenv_max')
	 params.tenv_max = 0.3;
end

if ~isfield(params, 'local_dc')
	params.local_dc = 0;
end

if ~isfield(params, 'normalize')
	 params.normalize = 1;
end

if ~isfield(params, 'show_or_preprocess')
	params.show_or_preprocess = 1;
end

if ~isfield(params, 'verbose')
	params.verbose = 1;
end

if ~isfield(params, 'valid_w_index')
	params.valid_w_index = NaN;
end

params.class = 'wavelets';

return;

%---------------------------------------------------------------------
% Making a list of gabor parameters
%---------------------------------------------------------------------
function gparams = getGaborParameters(params)

velangle_array = (0:params.veldivisions-1)/params.veldivisions*90;
sfmin_normalized = params.sfmin/params.sfmax;

if params.f_step_log
	freq_array = logspace(log10(sfmin_normalized), log10(1), params.fdivisions);
else
	freq_array = linspace(sfmin_normalized, 1, params.fdivisions);
end
dir_array = (0:params.dirdivisions-1)/params.dirdivisions*360;

switch params.phasemode
case 0  % spectral amplitudes
	pmarray = [0];
case 1  % linear sin and cos transform amplitudes
	pmarray = [1 2];
case 2  % half rectified sin and cos amplitudes
	pmarray = [3 4 5 6];
case 3  % 0+1
	pmarray = [0 1 2];
case 4  % 0+2
	pmarray = [0 3 4 5 6];
end

dirstart = 1;
if params.local_dc
	dirstart = 0; %% add local dc channels
end

waveletcount = 0;
gparams = zeros(8, 20000, 'single'); % prepare for some amount of memory for gparams

for vi=1:params.veldivisions
	velangle = velangle_array(vi);
	for fi = 1:params.fdivisions
		freq = freq_array(fi);
		sf = cos(velangle*pi/180)*params.sfmax*freq;
		tf = sin(velangle*pi/180)*params.tfmax*freq;
		senv = min([params.senv_max 1/sf*params.sf_gaussratio]);
		if tf ~= 0
			tenv = min([params.tenv_max 1/tf*params.tf_gaussratio]);
		else
			tenv = params.tenv_max;
		end
		numsps2 =floor((1-senv*params.std_step)/(params.std_step*senv)/2);
		numsps2 = max([numsps2 0]);
		centers = senv*params.std_step*(-numsps2:numsps2) + 0.5;
		[cx cy] = meshgrid(centers, centers);

		thisnumdirs = length(dir_array);
		if velangle == 0
			thisnumdirs = ceil(thisnumdirs/2);  % use only ~180 deg
		end
		if freq == 0
			thisnumdirs = 1;
		end
		for xyi = 1:length(cx(:))
			xcenter = cx(xyi);
			ycenter = cy(xyi);
			for diri = dirstart:length(dir_array)
				if diri
					dir = dir_array(diri); thissf = sf;
				else
					dir = 0; thissf = 0; % local dc channels
				end
				if  sf >= params.phasemode_sfmax
					waveletcount = waveletcount+1;
					thisgparam = [xcenter ycenter dir thissf tf senv tenv 0];
					gparams(:,waveletcount) = thisgparam;
				else
					for pmod = pmarray
						waveletcount = waveletcount+1;
						thisgparam = [xcenter ycenter dir thissf tf senv tenv pmod];
						gparams(:,waveletcount) = thisgparam;
					end
				end
			end
		end
	end
end

gparams = gparams(:,1:waveletcount);

