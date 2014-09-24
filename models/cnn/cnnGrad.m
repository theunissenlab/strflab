function [strf,g,savedStimFft]=nnGrad(strf,datIdx,savedStimFft)
%NNGRAD Evaluate gradient of error function for 2-layer network.
%

%
%       Description
%       G = NNGRAD(STRF, STIM, RESP) takes a generalized linear model data
%       structure STRF  together with a matrix STIM of input vectors and a
%       matrix RESP of responses, and evaluates the gradient G of the error
%       function with respect to the network weights. The error function
%       corresponds to the choice of output unit activation.
%
%
%       Usage:
%
%       G =NNGRAD(STRF, STIM, RESP)
%
%
%       Inputs:
%
%       STRF,   A STRF structure
%       STIM,   An nXm stimulus matrix (or vector)
%       RESP,   An nXtrials response matrix (or vector)
%   savedStimFft (optional) the Fourier transform of the stimulus, which is
%   sometimes saved and which makes Fourier domain calculations of the grad
%   much faster.
%
%       Outputs:
%
%       G, Value of the gradient of the error function for a given STRF,
%          STIM and RESP
%
%%
%      
%(Some code modified from NETLAB)


global globDat;

n_lag = length(strf.delays);
%  if n_lag > 1%(log(size(stim,1))/5) %Fourier-domain is faster here
if (~isfield(strf, 'freqDomain') & n_lag>30 ) | (isfield(strf, 'freqDomain') & strf.freqDomain)
    %disp('fd')
    if nargout == 1
        g=nnGradFourierDomain(strf,datIdx);
    else
        if exist('savedStimFft','var')
            if isnumeric(savedStimFft)
                g=nnGradFourierDomain(strf,datIdx,savedStimFft);
            else
                [g,savedStimFft]=nnGradFourierDomain(strf,datIdx);
            end
        else
            [g,savedStimFft]=nnGradFourierDomain(strf,datIdx);
        end
    end
else
    %disp('td')
    savedStimFft = 'Not used';
    g=nnGradTimeDomain(strf,datIdx);
end


[g, gprior] = nnGradBayes(strf, g);
