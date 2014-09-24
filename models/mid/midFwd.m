function [strf,resp_strf_dat,a_dat]=midFwd(strf,datIdx) %[resp_strf,a] = midFwd(strf, stim)
%
% function [[strf,resp_strf_dat,a_dat]=midFwd(strf,datIdx)
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
%	lin, midPak, midUnPak, midErr, linGrad
%
%
%(Some code modified from NETLAB)


% % Check arguments for consistency
% errstring = consist(strf, 'mid', x);
% if ~isempty(errstring);
%   error(errstring);
% end
% 
% a=s*strf.w1;

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

a = a/(max(a)-min(a));
% a = a + strf.b1;
resp_strf = a;
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




