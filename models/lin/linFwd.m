function [strf,resp_strf_dat,a_dat]=linFwd(strf,datIdx)
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
for ti=1:length(strf.delays)
  at = globDat.stim * strf.w1(:,ti);

  thisshift = strf.delays(ti);
  if thisshift>=0
    a(thisshift+1:end) = a(thisshift+1:end) + at(1:end-thisshift);
  else
    offset = mod(thisshift, samplesize);
    a(1:offset) = a(1:offset) + at(-thisshift+1:end);
  end

end
a = a + strf.b1;

switch strf.outputNL

  case 'linear'     % Linear outputs
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
nanmask = mod(strf.delays, size(globDat.stim,1)+1);
nanmask = nanmask(find(nanmask)); % no mask for delay 0
a(nanmask) = NaN;
resp_strf(nanmask) = NaN;


resp_strf_dat = resp_strf(datIdx);
a_dat = a(datIdx);

strf.internal.compFwd = 0;
strf.internal.prevResp = resp_strf;
strf.internal.prevLinResp = a;
strf.internal.dataHash = globDat.dataHash;
