function strf = midInit(nIn, delays, outfunc, numbins, binsize, startCond)
%GLMINIT	Create a generalized linear model.
%
%	Description:
%
%	STRF = MIDINIT(NIN, NOUT, OUTFUNC, DELAYS, NUMBINS, ) takes the number of inputs 
%	and outputs for a generalized linear model, together with a string 
%	OUTFUNC which specifies the output unit activation function, and 
%	returns a strf structure STRF. The weights are drawn from a zero mean,
%	isotropic Gaussian, with variance scaled by the fan-in of the output 
%	units. This makes use of the Matlab function RANDN and so the seed for 
%	the random weight initialization can be  set using RANDN('STATE', S)
%	where S is the seed value.
%	
%	Usage:
%
%	STRF = MIDINIT(NIN, OUTFUNC, DELAYS)
%
%
%	Inputs:
%
%	NIN, 	 Number of inputs in one time slice of the strf
%	OUTFUNC, String indicating the output function and regression type
%		 It can take the values:
%
%		
%               'nonparametric' (For use with a non-parametric
%               binned function )
%       
%	Outputs:
%
%	STRF, A STRF structure
%
%	The fields a glm strf are:
%
%	  type = 'glm'
%	  nin = number of inputs
%	  nwts = total number of weights and biases
%	  actfn = string describing the output unit activation function:
%	     
%             'nonparametric'
%	  w1 = first-layer weight matrix
%	  b1 = first-layer bias vector
%
%
%	See also
%	MIDPAK, MIDUNPAK, MIDFWD, MIDERR, MIDGRAD
%
%
%(Some code modified from NETLAB)

if nargin<2
  delays=[0];
end
if nargin<3
  outfunc='nonparametric';
end
if nargin<4
  numbins=100;
end
if nargin<5
  binsize=0.1;
end
if nargin<6
  startCond = 'const';
end
strf.type = 'mid';
strf.nIn = nIn;
strf.nWts = (nIn*length(delays) + 1);
strf.delays = delays;
strf.nbin = numbins;
strf.bin = binsize;

outtfns = {'nonparametric'};

if sum(strcmp(outfunc, outtfns)) == 0
  error('Undefined activation function. Exiting.');
else
  strf.outfn = outfunc;
end

% if nargin > 3
%   if isstruct(prior)
%     net.alpha = prior.alpha;
%     net.index = prior.index;
%   elseif size(prior) == [1 1]
%     net.alpha = prior;
%   else
%     error('prior must be a scalar or structure');
%   end
% end
  
% strf.w1 = randn(nIn, length(delays))/sqrt(nIn + 1);
% strf.b1 = randn(1)/sqrt(nIn + 1);
if strcmp(startCond, 'const')
	strf.w1 = ones(nIn, length(delays))/sqrt(nIn + 1);
	strf.b1 = ones(1)/sqrt(nIn + 1);
elseif strcmp(startCond, 'rand')
	strf.w1 = randn(nIn, length(delays))/sqrt(nIn + 1);
	strf.b1 = randn(1)/sqrt(nIn + 1);
end
	

% if nargin == 5
%   net.beta = beta;
% end

strf.internal.compFwd=1;
strf.internal.prevResp = [];
strf.internal.prevLinResp = [];
strf.internal.dataHash = NaN;
