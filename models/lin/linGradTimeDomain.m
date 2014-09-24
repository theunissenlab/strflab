function g=linGradTimeDomain(strf,datIdx)
%function [g] = linGrad(strf, datIdx)
%
% Evaluate gradient of error function for generalized linear model.
%
% Takes a generalized linear model data structure [strf]
% and evaluates the gradient [g] of the error
% function with respect to the network weights. The error function
% corresponds to the choice of output unit activation.  
% 
% 
% INPUT:	
%			[strf] = a linear STRF structure (see linInit for fields)
%			[datIdx] = a set of indices into the global [stim] matrix. 
%
% OUTPUT:
%			   [g] = a 1x(D*L) vector giving the gradient of the error function. 
%					 D=dimensions of stimulus, L=number of delays specified in strf.delays
%
%
%	SEE ALSO
%	linInit, linPak, linUnpak, linFwd, linErr, linGradFourierDomain, linGradTimeDomain
%
%(Some code modified from NETLAB)

global globDat;


%  resp_strf = linFwd(strf, 1:globDat.nSample);
[strf,resp_strf] = linFwd(strf, datIdx);

%  delout = resp_strf - globDat.resp;
if (size(resp_strf,1) == size(globDat.resp(datIdx),2)) & (size(resp_strf,2) == size(globDat.resp(datIdx),1))
    %This means the resps need to be transposed, that's all.
    resp_strf = resp_strf';
end


delout = resp_strf - globDat.resp(datIdx);

delays = strf.delays;
numdelays = length(delays);
delout(find(isnan(delout))) = 0;

%% assuming nout = 1
delout_array = repmat(delout(1)*0, [numdelays globDat.nSample]);

for ti=1:numdelays
	thisIdx = datIdx-delays(ti);
	validIdx = find(thisIdx>0 & thisIdx<=globDat.nSample);
	delout_array(ti,thisIdx(validIdx)) = delout(validIdx);
end

gw1 = (delout_array*globDat.stim)';
gb1 = nansum(delout);
g = gw1(:)';
if isfield(strf,'b1')
  g = [g gb1];
end
