function [stim, params] = preprocAnyFilters(rawStim, params)

optDeflt.class = 'preprocAnyFilters';
optDeflt.normalize = 1;
optRange.normalize = [0 1 2];
optDeflt.timeReverse = 0;
optRange.timeReverse = [0 1];
optDeflt.zscoreFilts = 0;
optRange.zscoreFilts = [0 1];
optDeflt.zeromean = 1;
optRange.zeromean = [0 1];
optDeflt.verbose = 1;
optRange.verbose = [0 1];
if nargin<2
  params=optDeflt;
else
  params=defaultOpt(params,optDeflt,optRange);
end
if nargin<1
  stim=optDeflt;
  return;
end

rawStim = single(rawStim);
stimxytsize = size(rawStim);
rawStim = reshape(rawStim, [prod(stimxytsize(1:end-1)) stimxytsize(end)]);

if params.zeromean
    if params.verbose, fprintf('[[zero mean stimuli]]\n'); end
	if isfield(params, 'zeromean_value')
		rawStim = rawStim - params.zeromean_value;
	else
		thismean = mean(rawStim(:));
		rawStim = rawStim - thismean;
		params.zeromean_value = thismean;
	end
end


params.filters = single(params.filters);
szFilts = size(params.filters);

if params.zscoreFilts
    % params.filters = norm_std_mean(reshape(params.filters, [prod(szFilts(1:end-1)) szFilts(end)]));
    params.filters = norm_std_mean(reshape(params.filters, [prod(szFilts(1:3)) size(params.filters,4)]));
end


params.filters = reshape(params.filters, [prod(szFilts(1:2)) szFilts(3) szFilts(4)]);

%% Time reverse filters if needed
if params.timeReverse
    params.filters = flipdim(params.filters, 2);
end


stim = zeros(stimxytsize(end), size(params.filters,3), 'single');

for ii = 1:size(params.filters,3)
    stim(:,ii) = dotdelay(squeeze(params.filters(:,:,ii)), rawStim);
    
    %% apply nonlinearity to filter output
    if isfield(params, 'nonlins')
        stim(:,ii) = feval(params.nonlins{ii}, stim(:,ii));
    end
    
end

if params.normalize==1
	if isfield(params, 'means') % Already preprocessed. Use the means and stds
		[stim] = norm_std_mean(stim, params.stds, params.means);
	else
		[stim, stds, means] = norm_std_mean(stim);
		params.means = means;
		params.stds = stds;
	end
elseif params.normalize==2
    if isfield(params, 'mins') % Already preprocessed. Use the mins and maxs
		[stim] = unit_scale(stim, params.mins, params.maxs);
	else
		[stim, mins, maxs] = unit_scale(stim);
		params.mins = mins;
		params.maxs = maxs;
	end
end

params.filters = reshape(params.filters, szFilts);