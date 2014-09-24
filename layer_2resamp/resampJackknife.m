function [strfArray,optArray,jackIdx]=resampJackknife(strf,datIdx,options)
%function [strfArray,optArray,jackIdx]=resampJackknife(strf,datIdx,options)
%
% Resample the data by holding a portion of the data out, and compute model
% STRF parameters from the remaining data.
%
% INPUT:
%      [strf] = model structure obtained via upper level *Init functions
%    [datIdx] = a vector containing indices of the samples to be use in the 
%               resampling fitting.
%   [options] = option structure containing fields:
%   .funcName = char array, name of this function.
%    .nResamp = # of re-fitting over the jackknifed samples.
%   .jackFrac = fraction of data to be held out.
% .useHoldOut = boolean, Default=0=Hold out data set is completely ignored.
%               1=allow to use hold out data set for early stopping or other 
%               indirect fitting functions. When this flag is 1, the .testFrac 
%               field is ignored.  (Note: using the hold out data for early
%               stopping goes against the principle of the Jackknife since no
%               data is actually held out of the fitting)
%   .testFrac = fraction of data to be use for early stopping or indirect 
%               fitting when the held out set is not allowed to be used.
%  .randomize = boolean, Default=1=randomize the order of the samples first
%               before splitting it n-fold. 0=Keep the original data order.
%   .randSeed = random seed used to initialize the state of rand.
%   .optimOpt = options required by the specified optimization algorithm.
%               Default uses trnGradDesc.m for training.
% OUTPUT:
% [strfArray] = structure array of structure containing [options.nResamp] 
%               fitted models.
%  [optArray] = structure array of options
%   [jackIdx] = Index of the hold out data set samples.
%
% SEE ALSO: resampBootstrap, resampCrossVal
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jul 2007)
% Bugfix for dealing with groups by Michael Oliver (Jan 2009)
% Change default to NOT use hold out set for early stopping (Oct 2009)
%
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='resampJackknife';
optRange.funcName={'resampJackknife'};

optDeflt.nResamp=10;
optRange.nResamp=[1,1e4];
optDeflt.jackFrac=0.10;
optRange.jackFrac=[0,1];
optDeflt.useHoldOut=0;
optRange.useHoldOut=[0,1];
optDeflt.testFrac=0.10;
optRange.testFrac=[0,1];
optDeflt.randomize=1;
optRange.randomize=[0,1];
optDeflt.randSeed=123456;
optRange.randSeed=[0,1e16];
optDeflt.optimOpt=feval('trnGradDesc');

if nargin<3
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end

if nargin<1
  strfArray=optDeflt;
  return;
end


if options.testFrac==0 & isfield(options.optimOpt, 'earlyStop') & options.optimOpt.earlyStop == 1
    warning(['You are using a fitting algorithm that uses an early stopping set, and have testFrac=0, therefore the held out' ...
        ' set of the jackknife will be used as the early stopping set!']);
    options.useHoldOut=1;
    pause(2);
elseif options.testFrac>0 & (~isfield(options.optimOpt, 'earlyStop') | options.optimOpt.earlyStop == 0) & options.useHoldOut==0;
    warning(['You are using a fitting algorithm that does not use an early stopping set, and have testFrac>0, therefore some' ...
        ' data is not being ignored from the training set']);
    pause(2);
end

% Get sample size & Check Group Index consistency
% --------------------

sampSiz=length(datIdx);
if ~isempty(globDat.groupIdx);
  groupIdx=globDat.groupIdx(datIdx);
  groupLen=length(groupIdx);
  if groupLen~=sampSiz
    error('resampJackknife >< group index must have same length as data index');
  end
  uniqGS=unique(groupIdx);
  nGroup=length(uniqGS);
else
  uniqGS=1:sampSiz;
  groupIdx=uniqGS;
  nGroup=sampSiz;
end


% Initial resampling computation & get model parameters
% --------------------
gHoldOutLen=round(options.jackFrac*nGroup);
gTrainLen=nGroup-gHoldOutLen;


% Compute models over Jackknife samples
% --------------------
rand('state',options.randSeed);


for ii=1:options.nResamp
  disp(sprintf('----- Jackknife Sample %d/%d',ii,options.nResamp));
  if options.randomize
    % rsIdx=randperm(nGroup);
    rsIdx = randsample(uniqGS,nGroup);
  else
    rsIdx=circshift(uniqGS,[0,(ii*gHoldOutLen)]);
  end
  
  holdOutIdx=findIdx(rsIdx(1:gHoldOutLen),groupIdx);
  trainGroup=rsIdx((gHoldOutLen+1):end);
  trainIdx=findIdx(trainGroup,groupIdx);
  
  if options.useHoldOut

    [strfArray(ii),optArray(ii)]=feval(options.optimOpt.funcName,strf, ...
      datIdx(trainIdx),options.optimOpt,datIdx(holdOutIdx));

  else  % don't useHoldOut
    testLen=round(options.testFrac*gTrainLen);
    testIdx=findIdx(trainGroup(1:testLen),groupIdx);
    fitIdx=findIdx(trainGroup((testLen+1):end),groupIdx);    
    
    [strfArray(ii),optArray(ii)]=feval(options.optimOpt.funcName,strf, ...
      datIdx(fitIdx),options.optimOpt,datIdx(testIdx));

  end  % options.useHoldOut
  
  jackIdx{ii}=holdOutIdx;
end


