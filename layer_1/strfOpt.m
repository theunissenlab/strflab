function [strf,options,varargout]=strfOpt(strf,datIdx,options,varargin)

%[strf,options,varargout]=strfOpt(strf,datIdx,options,varargin)
%
%	strfOpt is a helper function which facilitates the training of 
%	STRFs. It calls any of the middle-layer training functions to 
%	optimize the parameters, and passes a structure to govern the
%	optimization behaviour.
%
%
% INPUT:	
%	  [strf] = model structure obtained via upper level *Init functions
%   [datIdx] = a vector containing indices of the samples to be used in the 
%              fitting.
%  [options] = A structure which specifies and determines the behavior of the
%	       	   training function. (ex. options.funcName='resampBootstrap').
%		       Specific values for OPTIONS depend on which algorithm
%		       is used.  Type help trn* for that algorithm's
%		       options parameters.
%
% OUTPUT:
%
%    [strf] = structure containing model fit by the training algorithm
% [options] = option structure, with any additional fields added by fitting
%             algorithm
%
%
%(Some code modified from NETLAB)



[s{1:nargout}]=feval(options.funcName,strf,datIdx,options,varargin{:});
strf=s{1};

if nargout > 1
  options = s{2};

  % If there are additional arguments, extract them
  nextra = nargout - 2;
  if nextra > 0
    for ii = 1:nextra
      varargout{ii} = s{ii+2};
    end
  end
end
