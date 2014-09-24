function strf = cnnInit(filtDim,nHidden,frameSize,delays,outputNL,constrain,freqDom)
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
%	  filtDim = size and number of filters to fit i.e. [20 20 3]
%	  nhidden = number of hidden units between combined feature map(s) and output weights
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
if nargin<4
  delays=[0];
end
if nargin<5
  outputNL='linear';
end
if nargin<6
  constrain=0;
end
if nargin<7
  freqDom=0;
end


strf.type = 'cnn';
strf.filtDim = filtDim;
strf.nHidden = nHidden;
strf.constrain = constrain;
strf.frameSize = frameSize;
% strf.nWts = (nIn*length(delays) + 1)*nHidden + (nHidden + 1);

strf.delays = delays;

if constrain == 1
    nOutputMaps = filtDim(3);
    nFilt = 1;
elseif constrain == 2
    nOutputMaps = 1;
    nFilt = 1;
elseif constrain == 3
    nOutputMaps = filtDim(3);
    nFilt = length(delays);
else
    nOutputMaps = 1;
    nFilt = length(delays);
end

if nHidden > 0
    strf.nWts = prod([filtDim nFilt]) + (filtDim(3)*length(delays) + 1) + ((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1)+1)*nOutputMaps*nHidden + (nHidden+1);
    strf.hlw = randn((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1)*nOutputMaps,nHidden)/sqrt((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1) + 1);
    strf.hlb = randn(1, nHidden)/sqrt((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1) + 1);
    strf.outw = randn(nHidden, 1)/sqrt(nHidden + 1);
    strf.outb = randn(1)/sqrt(nHidden + 1);
else
    strf.nWts = prod([filtDim nFilt]) + (filtDim(3)*length(delays) + 1) + ((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1)+1);
    strf.hlw = [];
    strf.hlb = [];
    strf.outw = randn((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1)*nOutputMaps, 1)/sqrt((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1) + 1);
    strf.outb = randn(1)/sqrt((frameSize(1)-filtDim(1)+1)*(frameSize(2)-filtDim(2)+1) + 1);
end

% strf.w1 = randn(nIn,nHidden,length(delays))/sqrt(nIn*length(delays) + 1);
strf.filts = randn([filtDim nFilt])/sqrt(prod([filtDim(1:2) nFilt])+1);
strf.fmw = randn([filtDim(3) length(delays)])/sqrt(filtDim(3)*nFilt+1);
strf.fmb = randn(1, nOutputMaps)/sqrt(prod(frameSize)+1);
% strf.b1 = randn(1, nHidden)/sqrt(nIn + 1);
% strf.w2 = randn(nHidden, 1)/sqrt(nHidden + 1);
% strf.b2 = randn(1)/sqrt(nHidden + 1);

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


