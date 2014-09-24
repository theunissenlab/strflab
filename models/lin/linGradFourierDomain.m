function [gOut,savedStimFft] = linGradFourierDomain(strf, stim, resp,savedStimFft)
%function [g,savedStimFft] = linGrad(strf, savedStimFft)
%
% Evaluate gradient of error function for generalized linear model using Fourier
% Domain representation.
%
% Takes a generalized linear model data structure [strf]
% and evaluates the gradient [g] of the error
% function with respect to the network weights. The error function
% corresponds to the choice of output unit activation.  This function may
% be called by linGrad by setting strf.fourierdoman = 1.  Users need not call this
% function directly.
% 
% 
%
% INPUT:	
%			[strf] = a linear STRF structure (see linInit for fields)
%   [savedStimFft] = (optional) the Fourier transform of the stimulus, which is
%   				 sometimes saved and which makes Fourier domain calculations of the grad
%   				 much faster.
%
% OUTPUT:
%			   [g] = a 1x(D*L) vector giving the gradient of the error function. 
%					 D=dimensions of stimulus, L=number of delays specified in strf.delays
%   [savedStimFft] = the Fourier transform of the stimulus.
%					 
%
%
%	SEE ALSO
%	linInit, linPak, linUnpak, linFwd, linErr, linGrad
%
%(Some code modified from NETLAB)



resp_strf = linFwd(strf, stim);

delout_all = resp_strf - resp;

%numdelays = length(strf.delays);
delays = strf.delays;
maxlag = max(delays);
% gw1 = zeros(size(stim,2), strf.nout, numdelays);
% for ti=1:numdelays
%         sdelout = shift(delout, [-strf.delays(ti) 0]);
%         nonnans = find(~isnan(sdelout(:,1)));
%         gw1(:,:,ti) = stim(nonnans,:)'*sdelout(nonnans,:);
% end

gOut = [];

    delout = delout_all(:,1);
    if nargout > 1
        if ~exist('savedStimFft','var')
            [g2,savedStimFft] = stim_error_corr(stim,delout,maxlag);
        else
            [g2,savedStimFft] = stim_error_corr(stim,delout,maxlag,savedStimFft);
        end
    else
        if ~exist('savedStimFft','var')
            [g2] = stim_error_corr(stim,delout,maxlag);
        else
            [g2] = stim_error_corr(stim,delout,maxlag,savedStimFft);
        end
    end

    S = size(stim,2);
    g = zeros(1, length(delays) * S);
    for jj = 1:length(delays)
        g((1:S) + (jj-1)*S) = g2(maxlag+1 - delays(jj),:);
    end
    gOut = [gOut g];

gb1 = nansum(delout);
if isfield(strf,'b1')
  gOut = [gOut gb1];
end
