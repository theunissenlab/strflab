function strf = nnInit(nIn, nHidden,delays,outputNL,freqDom)
%function strf = nnInit(nIn, nHidden,delays,outputNL,freqDom)
%
%NNINIT	Create a 2-layer feedforward network.
%
%	Description
%	STRF = NNINIT(NIN, NHIDDEN, FUNC, DELAYS) takes the number of inputs,
%	hidden units and output units for a 2-layer feed-forward network,
%	together with a string FUNC which specifies the output unit
%	activation function, and returns a data structure STRF. The weights
%	are drawn from a zero mean, unit variance isotropic Gaussian, with
%	varianced scaled by the fan-in of the hidden or output units as
%	appropriate. This makes use of the Matlab function RANDN and so the
%	seed for the random weight initialization can be  set using
%	RANDN('STATE', S) where S is the seed value.  The hidden units use
%	the TANH activation function.
%
%	The fields in NET are
%	  type = 'mlp'
%	  nin = number of inputs
%	  nhidden = number of hidden units
%	  nwts = total number of weights and biases
%	  actfn = string describing the output unit activation function:
%	      'linear'
%	      'logistic
%	      'softmax'
%	  w1 = first-layer weight matrix
%	  b1 = first-layer bias vector
%	  w2 = second-layer weight matrix
%	  b2 = second-layer bias scalar
%	  delays = vector of time delays
%	 Here W1 has dimensions NIN times NHIDDEN, B1 has dimensions 1 times
%	NHIDDEN, W2 has dimensions NHIDDEN times 1
%
%Some code taken from Netlab
%

% Set default option values
% --------------------
if nargin<3
  delays=[0];
end
if nargin<4
  outputNL='linear';
end
if nargin<5
  freqDom=0;
end

strf.type = 'nn';
strf.nIn = nIn;
strf.nHidden = nHidden;
strf.nWts = (nIn*length(delays) + 1)*nHidden + (nHidden + 1);
strf.delays = delays;

strf.w1 = randn(nIn,nHidden,length(delays))/sqrt(nIn*length(delays) + 1);
strf.b1 = randn(1, nHidden)/sqrt(nIn + 1);
strf.w2 = randn(nHidden, 1)/sqrt(nHidden + 1);
strf.b2 = randn(1)/sqrt(nHidden + 1);

nlSet={'linear', 'logistic', 'exponential', 'softmax'};
if ismember(outputNL,nlSet)
  strf.outputNL = outputNL;
else
  error('nnInit >< Unknown Output Nonlinearity!');
end
strf.freqDomain = freqDom;

strf.internal.compFwd=1;
strf.internal.prevResp = [];
strf.internal.prevLinResp = [];
strf.internal.prevHidResp = [];
strf.internal.dataHash = NaN;


