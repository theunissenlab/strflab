function [strf,g,savedStimFft]=linGrad(strf,datIdx,savedStimFft)
%function [strf, g,savedStimFft] = linGrad(strf, datIdx, savedStimFft)
%
% Evaluate gradient of error function for generalized linear model.
%
% Takes a generalized linear model data structure [strf]
% and evaluates the gradient [g] of the error
% function with respect to the network weights. The error function
% corresponds to the choice of output unit activation.  If [savedStimFft]
% is included as an input, the gradient will be evaluated using a fourier domain 
% technique.  Alternatively, this method may be chosen by setting strf.fourierdoman = 1.
% 
% 
%
% INPUT:	
%			[strf] = a linear STRF structure (see linInit for fields)
%			[datIdx] = a set of indices into the global [stim] matrix. 
%   [savedStimFft] = (optional) the Fourier transform of the stimulus, which is
%   				 sometimes saved and which makes Fourier domain calculations of the grad
%   				 much faster.
%
% OUTPUT:
%			[strf] = unmodified STRF structure.
%			   [g] = a 1x(D*L) vector giving the gradient of the error function. 
%					 D=dimensions of stimulus, L=number of delays specified in strf.delays
%   [savedStimFft] = (optional) if the fourier-based gradient has been selected, will contain
%					 the Fourier transform of the stimulus.
%
%
%	SEE ALSO
%	linInit, linPak, linUnpak, linFwd, linErr, linGradFourierDomain, linGradTimeDomain
%
%(Some code modified from NETLAB)

global globDat;

n_lag = length(strf.delays);
%  if n_lag > 1%(log(size(stim,1))/5) %Fourier-domain is faster here
if (~isfield(strf, 'freqDomain') & n_lag>30 ) | (isfield(strf, 'freqDomain') & strf.freqDomain)
    %disp('fd')
    if nargout == 1
        g=linGradFourierDomain(strf,datIdx);
    else
        if exist('savedStimFft','var')
            if isnumeric(savedStimFft)
                g=linGradFourierDomain(strf,datIdx,savedStimFft);
            else
                [g,savedStimFft]=linGradFourierDomain(strf,datIdx);
            end
        else
            [g,savedStimFft]=linGradFourierDomain(strf,datIdx);
        end
    end
else
    %disp('td')
    savedStimFft = 'Not used';
    g=linGradTimeDomain(strf,datIdx);
end

