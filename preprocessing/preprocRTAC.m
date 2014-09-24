function [stim, params] = preprocRTAC(rawStim, params)

% function [stim, params] = preprocRTAC(rawStim, params);
%
% A script for preprocessing of stimuli in order to perform RTA and/or RTC
%
% INPUT:
%           [rawStim] = A X-by-Y-by-T matrix containing stimuli (movie)
%            [params] = structure that contains parameters for preprocessing
%            .covtime = determines whether to compute covariance across time (default: 0)
%               .RTAC = vector determined whether to preprocess for RTA and/or RTC
%                       i.e. [1 0] is just RTA, [1 1] is both RTA and RTC (default: [1 0])
%          .covdelays = number of frames back to compute RTC for when calculating covariance over time
%                       (default: 0)
%           .locality = parameter that determines the spatial extent of pixel interactions
%                       (default: 1)
%          .normalize = A flag to normalize mean and stds for each channel (default: 1)
%
% OUTPUT:
%              [stim] = Preprocessed stimuli that can be used for STRF fitting.
%                       NxD matrix, N=sample size, D=dimensions.
%            [params] = structure that contains parameters for preprocessing, with additional fields:
%
%
% EXAMPLE:
%  params = preprocRTAC;
%    returns the default set of parameters.
%
%  [stim, params] = preprocRTAC(rawStim, params)
%    returns preprocessed stimuli and parameters
%
% SEE ALSO: preprocWavelets3d, preprocSpectra
% ====================
% Edit February 2010 to normalize all channels

optDeflt.class = 'preprocRTAC';
optDeflt.covtime = 0;
optRange.covtime = [0 1];
optDeflt.RTAC = [1 0];
optRange.RTAC = [0 1];
optDeflt.covdelays = 0;
optRange.covdelays = [0 Inf];
optDeflt.locality = 1;
optRange.locality = [0 Inf];
optDeflt.normalize = 1;
optRange.normalize = [0 1];
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
params.covsz = size(rawStim,1)*size(rawStim,2);

if ~params.covtime
    params.covdelays = 0;
end

if params.covtime & params.covdelays == 0
    error('You should set the covdelays to greater than zero to calculate covariance across time');
end

if params.RTAC(1) & params.RTAC(2)
    [A, idxr, idxc] = getSubMatrix([size(rawStim,1) size(rawStim,2) (params.covdelays + 1)], params.locality*2 +1, [1 1 params.covtime]);
    params.idxr = idxr; params.idxc = idxc;
    rawStim = reshape(rawStim, size(rawStim,1)*size(rawStim,2), size(rawStim,3));
    
    if ~params.covtime
        stim = zeros([size(rawStim,2) size(idxr,1)+params.covsz], 'single');
        stim(:,1:params.covsz) = rawStim';
        %% Calculate the mean to subtract off
        ErawStim = mean(rawStim,2);
        rawStim = bsxfun(@minus,rawStim,ErawStim);
        stim(:,params.covsz+1:end) = (rawStim(idxr,:).*rawStim(idxc,:))';
    
    else
        stim = zeros([size(rawStim,2) size(idxr,1)+((params.covdelays + 1)*params.covsz)], 'single');
        for ii = 1:(params.covdelays + 1)
            stim(:,(ii-1)*params.covsz+1:ii*params.covsz) = [zeros([ii-1 params.covsz], 'single'); rawStim(:,1:end-ii+1)'];
        end
        
        %% Calculate the mean to subtract off
        ErawStim = mean(rawStim,2);
        rawStim = [zeros(size(rawStim,1), (params.covdelays + 1)-1) rawStim];
        rawStim = bsxfun(@minus,rawStim,ErawStim);
        
        for ii = 1:1000:size(idxr,1)
            iiend = min([size(idxr,1) ii+999]);
            %% Separate row and column indexes into pixel and time indicies
            [rr,tr] = ind2sub([size(rawStim,1) (params.covdelays + 1)], idxr(ii:iiend));
            [cc,tc] = ind2sub([size(rawStim,1) (params.covdelays + 1)], idxc(ii:iiend));
            
            %% Create pixel and time indices across all samples
            rridx = repmat(rr, [1 size(rawStim,2)-(params.covdelays + 1)+1]);
            ccidx = repmat(cc, [1 size(rawStim,2)-(params.covdelays + 1)+1]);
            
            tridx = bsxfun(@minus, repmat(uint32([(params.covdelays + 1):size(rawStim,2)]), [size(rr,1) 1]),tr)+1;
            tcidx = bsxfun(@minus, repmat(uint32([(params.covdelays + 1):size(rawStim,2)]), [size(cc,1) 1]),tc)+1;
            
            %% Get linear index into stimuli from pixel and time indicies
            indr = sub2ind(size(rawStim), rridx(:), tridx(:));
            indc = sub2ind(size(rawStim), ccidx(:), tcidx(:));
            
            stim(:,(params.covdelays + 1)*params.covsz+ii:(params.covdelays + 1)*params.covsz+iiend) = reshape(rawStim(indr).*rawStim(indc), [iiend-ii+1 size(rawStim,2)-(params.covdelays + 1)+1])';
            fprintf('preprocRTAC percent complete: %.1f\n', 100*iiend/size(idxr,1));
        end
    end
elseif params.RTAC(1) & ~params.RTAC(2)
    rawStim = reshape(rawStim, size(rawStim,1)*size(rawStim,2), size(rawStim,3));
    
    % stim = rawStim';
    stim = zeros([size(rawStim,2) ((params.covdelays + 1)*params.covsz)], 'single');
    for ii = 1:(params.covdelays + 1)
        stim(:,(ii-1)*params.covsz+1:ii*params.covsz) = [zeros([ii-1 params.covsz], 'single'); rawStim(:,1:end-ii+1)'];
    end
    
elseif ~params.RTAC(1) & params.RTAC(2)
    [A, idxr, idxc] = getSubMatrix([size(rawStim,1) size(rawStim,2) (params.covdelays + 1)], params.locality*2 +1, [1 1 params.covtime]);
    rawStim = reshape(rawStim, size(rawStim,1)*size(rawStim,2), size(rawStim,3));
    
    params.idxr = idxr; params.idxc = idxc;
    
    if ~params.covtime        
        %% Calculate the mean to subtract off
        ErawStim = mean(rawStim,2);
        rawStim = bsxfun(@minus,rawStim,ErawStim);
        stim= (rawStim(idxr,:).*rawStim(idxc,:))';
    else
        %% Calculate the mean to subtract off
        ErawStim = mean(rawStim,2);
        stim = zeros([size(rawStim,2) size(idxr,1)], 'single');
        rawStim = [zeros(size(rawStim,1), (params.covdelays + 1)-1) rawStim];
        rawStim = bsxfun(@minus,rawStim,ErawStim);  
        
        for ii = 1:1000:size(idxr,1)
            iiend = min([size(idxr,1) ii+999]);
            %% Separate row and column indexes into pixel and time indicies
            [rr,tr] = ind2sub([size(rawStim,1) (params.covdelays + 1)], idxr(ii:iiend));
            [cc,tc] = ind2sub([size(rawStim,1) (params.covdelays + 1)], idxc(ii:iiend));
            
            %% Create pixel and time indices across all samples
            rridx = repmat(rr, [1 size(rawStim,2)-(params.covdelays + 1)+1]);
            ccidx = repmat(cc, [1 size(rawStim,2)-(params.covdelays + 1)+1]);
            
            tridx = bsxfun(@minus, repmat(uint32([(params.covdelays + 1):size(rawStim,2)]), [size(rr,1) 1]),tr)+1;
            tcidx = bsxfun(@minus, repmat(uint32([(params.covdelays + 1):size(rawStim,2)]), [size(cc,1) 1]),tc)+1;
            
            %% Get linear index into stimuli from pixel and time indicies
            indr = sub2ind(size(rawStim), rridx(:), tridx(:));
            indc = sub2ind(size(rawStim), ccidx(:), tcidx(:));

            %% Calculate 2nd order terms
            stim(:, ii:iiend) = reshape(rawStim(indr).*rawStim(indc), [iiend-ii+1 size(rawStim,2)-(params.covdelays + 1)+1])';
            fprintf('preprocRTAC percent complete: %.1f\n', 100*iiend/size(idxr,1));
        end
    end
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