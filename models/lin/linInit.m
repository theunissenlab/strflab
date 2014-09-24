function strf=linInit(nIn,delays,outputNL,freqDom)
%function [strf] = linInit(nin, delays, outputNL, freqDom)
%
%Create a generalized linear model.
%
% Takes the number of inputs and outputs for a generalized linear model, together 
% with a string [outfunc] which specifies the output unit activation function, and 
% returns a strf structure [strf]. The weights are drawn from a zero mean,
% isotropic Gaussian, with variance scaled by the fan-in of the output 
% units. This makes use of the Matlab function RANDN and so the seed for 
% the random weight initialization can be  set using RANDN('STATE', S)
% where S is the seed value.
%
% INPUT:
% 		[nin] = number of inputs in one time slice of the strf
%    [delays] = row vector of integers, describing latencies for the strf.  
%		        (For example [0 1 2] is 3 time lags including no lag)
%	[outputNL] = string indicating the output function and regression type:
%				'linear'     (For use with squared error)
%				'logistic'   (For use with logistic regression)
%				'exponential'(For use with poisson regression)
%   [freqDom] = toggles between options for convolution of stimulus with model weight.  For models with long delays
%				(e.g., auditory data), try freqDom = 1.  For models with shorter delays,
%				freqDom = 0 will usually produce faster results.
%
% OUTPUT:
%	[strf] = a STRF structure containing fields:
%	   .type  = 'glm'
%	  	.nin  = number of inputs in one time slice of the strf (see above)
%  	 .delays  = delay vector input (see above)
%	   .nWts  = total number of weights and biases
%	  .actfn  = string describing the output unit activation function: 'linear',
%			    'logistic', or 'exponential'
% .freqDomain = 1 for convolution in frequency domain.
%	   	  .w1 = nin x 1 vector of weights
%	      .b1 = scalar bias term
%	
%
%
%	SEE ALSO
%	linPak, linUnpak, linFwd, linErr, linGrad
%
%(Some code modified from NETLAB)


% Set default option values
% --------------------
if nargin<2
  delays=[0];
end
if nargin<3
  outputNL='linear';
end
if nargin<4
  freqDom=0;
end

strf.type = 'lin';
strf.nIn = nIn;
strf.nWts = (nIn*length(delays) + 1);
strf.delays = delays;

%strf.w1 = randn(nIn, 1, length(delays))/sqrt(nIn + 1);
%strf.b1 = randn(1,1)/sqrt(nIn + 1);
strf.w1=zeros(nIn,length(delays));
strf.b1=0;

nlSet={'linear','logistic','softmax','exponential'}; 
if ismember(outputNL,nlSet)
  strf.outputNL = outputNL;
else
  error('linInit >< Unknown Output Nonlinearity!');
end
strf.freqDomain = freqDom;

strf.internal.compFwd=1;
strf.internal.prevResp = [];
strf.internal.prevLinResp = [];
strf.internal.dataHash = NaN;
