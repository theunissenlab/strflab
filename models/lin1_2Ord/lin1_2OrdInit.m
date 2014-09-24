function strf=lin1_2OrdInit(nIn,delays1,delays2,outputNL,interactIdx,freqDom)
%function [strf] = lin2ndOrdInit(nin, delays, outputNL, freqDom)
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
  delays1=[0];
  delays2=[0];
end
if nargin<4
  outputNL='linear';
end
if nargin<5
    interactIdx.idxr = [];
    interactIdx.idxc = [];
end
if nargin<6
  freqDom=0;
end
global globDat;

strf.type = 'lin1_2Ord';
strf.nIn = nIn;
strf.nWts = (size(globDat.stim,2)*length(delays1) + nIn*length(delays2) + 1);
strf.delays1 = delays1;
strf.delays2 = delays2;

%strf.w1 = randn(nIn, 1, length(delays))/sqrt(nIn + 1);
%strf.b1 = randn(1,1)/sqrt(nIn + 1);


nlSet={'linear','logistic','softmax','exponential','huber','logcosh'}; 
if ismember(outputNL,nlSet)
  strf.outputNL = outputNL;
else
  error('linInit >< Unknown Output Nonlinearity!');
end
strf.freqDomain = freqDom;
strf.interactIdx = interactIdx;
strf.normalize = 1;
strf.covdelays=max([interactIdx.idxrd; interactIdx.idxcd])-1;
strf.chunksize = 1000;
strf.internal.compFwd=1;
strf.internal.prevResp = [];
strf.internal.prevLinResp = [];
strf.internal.dataHash = NaN;


strf.w1=zeros(size(globDat.stim,2),length(delays1));
strf.w2=zeros(nIn,length(delays2));
strf.b1=0;



if size(globDat.stim,1)-strf.covdelays ~= size(globDat.resp,1)
    globDat.stim = [zeros(strf.covdelays, size(globDat.stim,2)); globDat.stim];
    
    % Compute hash for the data set 
    % --------------------
    respHash = 100*abs(nanmean(double(globDat.resp(:))) + nanmean(double(globDat.resp(1:11:end))));
    stimHash = 100*abs(nanmean(double(globDat.stim(1:109:end))));
    magdif = log10((respHash+0.00012)/(stimHash+0.00011));
    dataHash = respHash + stimHash * 10^magdif;
    globDat.dataHash = dataHash;
end




