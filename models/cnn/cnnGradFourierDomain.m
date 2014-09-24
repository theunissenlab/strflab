function [gOut,savedStimFft] = nnGradFourierDomain(strf, stim, resp,savedStimFft)

%NNGRADFOURIERDOMAIN Evaluate gradient of error function for neural network.
%		     (Using fourier transform method)
%
%	Description
%	G = NNGRADFOURIERDOMAIN(STRF, STIM, RESP) takes a neural network strf 
%	model data structure STRF  together with a matrix STIM of input vectors 
%	and a matrix RESP of responses, and evaluates the gradient G of the error
%	function with respect to the network weights, using the fourier transform
%
%
%	Usage:
%	
%	G = NNGRADFOURIERDOMAIN(STRF, STIM, RESP)
%
%
%	Inputs:
%		
%	STRF,	A STRF structure
%	STIM,	An nXm stimulus matrix (or vector)
%	RESP,	An nXtrials response matrix (or vector)
%
%
%	Outputs:
%
%	G, Value of the gradient of the error function for a given STRF,
%	   STIM and RESP
%
%
%
%
%(Some code modified from NETLAB)


[resp_strf, z] = nnFwd(strf, stim);

delout = resp_strf - resp;

gw2 = z'*delout;
gb2 = sum(delout, 1);

% Now do the backpropagation.
delhid = delout*strf.w2';

delhid = delhid.*(1.0 - z.*z);

gb1 = sum(delhid, 1);

delays = strf.delays;
maxlag = max(delays);

gOut = [];

for ii=1:strf.nhidden

    if nargout > 1
        if ~exist('savedStimFft','var')
            [g2,savedStimFft] = stim_error_corr(stim,delhid(:,ii),maxlag);
        else
            [g2,savedStimFft] = stim_error_corr(stim,delhid(:,ii),maxlag,savedStimFft);
        end
    else
        if ~exist('savedStimFft','var')
            [g2] = stim_error_corr(stim,delhid(:,ii),maxlag);
        else
            [g2] = stim_error_corr(stim,delhid(:,ii),maxlag,savedStimFft);
        end
    end

    S = size(stim,2);
    g = zeros(1, length(delays) * S);
    for jj = 1:length(delays)
        g((1:S) + (jj-1)*S) = g2(maxlag+1 - delays(jj),:);
    end
    gOut = [gOut g];

end

gOut = [gOut gb1];

gOut = [gOut, gw2(:)'];

gOut = [gOut, gb2];
