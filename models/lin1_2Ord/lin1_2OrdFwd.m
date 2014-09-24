function [strf,resp_strf_dat,a_dat]=lin2ndOrdFwd(strf,datIdx)
% function [strf,resp_strf_dat,a_dat] = linFwdFaster(strf, datIdx)
%  
% Forward propagation through generalized linear model.
%	
% Takes a generalized linear model data structure
% [strf] together with a stimulus matrix [stim] and forward 
% propagates the inputs through the network to generate a vector 
% [resp_strf] of output vectors.  Can also return [a], the vector of 
%  model responses generated before application of the output nonlinearity.
%
% INPUTS:
%
%	[strf]   = a strf stucture
%   [datIdx] = a set of indices into the global [stim] matrix.  Picks out samples to be propagated through model.
%
%
% OUTPUTS:
%			 [strf] = unmodified strf structure
%	[resp_strf_dat] = Nx1 vector of model outputs elicited by [strf] structure by stimuli identified in datIdx.
%			[a_dat] = Nx1 vector of response before output nonlinearity.
%
%
%	SEE ALSO
%	lin, linPak, linUnPak, linErr, linGrad
%
%
%(Some code modified from NETLAB)

global globDat;

samplesize = globDat.nSample;

if strf.internal.compFwd == 0 & samplesize == length(strf.internal.prevResp) & strf.internal.dataHash == globDat.dataHash
	resp_strf_dat = strf.internal.prevResp(datIdx);
	a_dat = strf.internal.prevLinResp(datIdx);
	return
end


a = zeros(samplesize, 1);

if ~isfield(strf, 'means1') && strf.normalize
    strf.means1 = zeros(1,size(globDat.stim,2));
    strf.stds1 = zeros(1,size(globDat.stim,2));
end

if ~isfield(strf, 'means2') && strf.normalize
    strf.means2 = zeros(size(strf.interactIdx.idxr))';
    strf.stds2 = zeros(size(strf.interactIdx.idxr))';
end

%% First order terms

for ti=1:length(strf.delays1)
    at = zeros(samplesize, 1);
    for ii = 1:strf.chunksize:size(globDat.stim,2)
        iiend = min([size(globDat.stim,2) ii+strf.chunksize-1]);
        stim = zeros(samplesize, iiend-ii+1);
        if strf.normalize
            if strf.stds1(ii) ~= 0
                [stim] = norm_std_mean(globDat.stim(strf.covdelays+1:end,ii:iiend), strf.stds1(ii:iiend), strf.means1(ii:iiend));
            else
                [stim, strf.stds1(ii:iiend), strf.means1(ii:iiend)] = norm_std_mean(globDat.stim(strf.covdelays+1:end,ii:iiend));
            end
        else
            [stim] = globDat.stim(strf.covdelays+1:end,ii:iiend);
        end
        at = at + (stim * strf.w1(ii:iiend,ti));
    end

  thisshift = strf.delays1(ti);
  if thisshift>=0
    a(thisshift+1:end) = a(thisshift+1:end) + at(1:end-thisshift);
  else
    offset = mod(thisshift, samplesize);
    a(1:offset) = a(1:offset) + at(-thisshift+1:end);
  end

end

%% Second order terms

for ti=1:length(strf.delays2)
    
  at = zeros(samplesize, 1);
  for ii = 1:strf.chunksize:size(strf.interactIdx.idxr,1)
      iiend = min([size(strf.interactIdx.idxr,1) ii+strf.chunksize-1]);
      stim = zeros(samplesize, iiend-ii+1);
      if strf.normalize
          if strf.stds2(ii) ~= 0
              for rr = 1:strf.covdelays + 1
                  for cc=1:strf.covdelays + 1
                      idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                      if ~isempty(idxf)
                          [stim(:,(idxf-ii+1))] = norm_std_mean(globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                          globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf)), strf.stds2(idxf), strf.means2(idxf));
                      end
                  end
              end
              
              
          else ...
              for rr = 1:strf.covdelays + 1
                  for cc=1:strf.covdelays + 1
                      idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                      if ~isempty(idxf)
                          [stim(:,(idxf-ii+1)), strf.stds2(idxf), strf.means2(idxf)] = norm_std_mean(globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                            globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf)));
                      end
                  end
              end
              
              
          end
      else ...
          for rr = 1:strf.covdelays + 1
              for cc=1:strf.covdelays + 1
                  idxf = intersect(find(strf.interactIdx.idxrd(ii:iiend)==rr),find(strf.interactIdx.idxcd(ii:iiend)==cc))+ii-1;
                  if ~isempty(idxf)
                      [stim(:,(idxf-ii+1))] = globDat.stim((strf.covdelays + 2 - rr):(end-rr+1),strf.interactIdx.idxr(idxf)).*...
                        globDat.stim((strf.covdelays + 2 - cc):(end-cc+1),strf.interactIdx.idxc(idxf));
                  end
              end
          end
          
      end
      at = at + (stim * strf.w2(ii:iiend,ti));
      fprintf('.')
  end

  thisshift = strf.delays2(ti);
  if thisshift>=0
    a(thisshift+1:end) = a(thisshift+1:end) + at(1:end-thisshift);
  else
    offset = mod(thisshift, samplesize);
    a(1:offset) = a(1:offset) + at(-thisshift+1:end);
  end

end

fprintf(' strfFwd step done.\n')

if isfield(strf,'b1')
    a = a + strf.b1;
end


switch strf.outputNL

  case 'linear'     % Linear outputs
    resp_strf = a;
    
  case 'huber'     % Linear outputs
    resp_strf = a;

  case 'logcosh'     % Linear outputs
    resp_strf = a;

  case 'logistic'   % Logistic outputs
    % Prevent overflow and underflow: use same bounds as glmerr
    % Ensure that log(1-y) is computable: need exp(a) > eps
    maxcut = -log(eps);
    % Ensure that log(y) is computable
    mincut = -log(1/realmin - 1);
    a = min(a, maxcut);
    a = max(a, mincut);
    resp_strf = 1./(1 + exp(-a));

  case 'softmax'        % Softmax outputs
    nout = size(a,2);
    % Prevent overflow and underflow: use same bounds as glmerr
    % Ensure that sum(exp(a), 2) does not overflow
    maxcut = log(realmax) - log(nout);
    % Ensure that exp(a) > 0
    mincut = log(realmin);
    a = min(a, maxcut);
    a = max(a, mincut);
    temp = exp(a);
    resp_strf = temp./(sum(temp, 2)*ones(1,nout));
    % Ensure that log(y) is computable
    resp_strf(resp_stf<realmin) = realmin;

  case 'exponential'
    resp_strf=exp(a);
  otherwise
    error(['Unknown activation function ', strf.outputNL]);
end

% mask for nonvalid frames
nanmask = mod(strf.delays1, size(globDat.stim(strf.covdelays+1:end,:),1)+1);
nanmask = nanmask(find(nanmask)); % no mask for delay 0
a(nanmask) = NaN;
resp_strf(nanmask) = NaN;


resp_strf_dat = resp_strf(datIdx);
a_dat = a(datIdx);

strf.internal.compFwd = 0;
strf.internal.prevResp = resp_strf;
strf.internal.prevLinResp = a;
strf.internal.dataHash = globDat.dataHash;
