function [strf,options]=trnBackslash(strf,datIdx,options,varargin)

%function [strf,options]=trnSCG(strf,datIdx,options)
%
% Scaled conjugate gradient optimization of STRF
%
% INPUT:
%            [strf] = model structure obtained via upper level *Init functions
%          [datIdx] = a vector containing indices of the samples to be used in the 
%                     fitting.
%         [options] = option structure containing fields:
%         .funcName = char array, name of this function.
%          .display = Number of iterations between displaying of errors as we train.
%                     Positive integer means display in stdout.
%                     Negative integer means display in a figure window.
%                     0 means do not display errors (default)
%          .maxIter = Maximum number of training iterations (default: 5000)
%     .minErrChange = The (positive real) size of error change which determines convergence
%	                  (default = eps(class(stim)))
%      .minStepSize = The (positive real) size of step size which determines convergence
%	                  (default = eps(class(stim)))
%
%
% OUTPUT:
%            [strf] = structure containing model fit by trnGradDesc.
%         [options] = option structure with additional fields:
%      .diagnostics = contains information on the fitting process.
% .diagnostics.errs = this has dimensions 1 x N where this has the error on the
%                     training set, over the course of the fitting iterations.
%
% EXAMPLE:
% options = trnSCG
%   returns default options
%
%
% SEE ALSO: trnGradDesc, trnDirectFit
%
%(Some code modified from NETLAB)
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='trnBackslash';
optRange.funcName={'trnBackslash'};


if nargin<3
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end
if nargin<1
  strf=optDeflt;
  return;
end

x = [globDat.stim(datIdx,:) ones([size(globDat.stim(datIdx,:),1) 1])] \ globDat.resp(datIdx);

strf=strfUnpak(strf,x');
