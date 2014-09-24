function [strfArray,optArray,bootIdx]=resampBootstrap(strf,datIdx,options)
%function [strfArray,optArray,bootIdx]=resampBootstrap(strf,datIdx,options)
%
% Resample the data with (or without) replacement and compute model STRF parameters
% for each set of resampled data.
%
% INPUT:
%      [strf] = model structure obtained via upper level *Init functions
%    [datIdx] = a vector containing indices of the samples to be used in the 
%               resampling fitting.
%   [options] = option structure containing fields:
%   .funcName = char array, name of this function.
%    .nResamp = # of re-fitting over the bootstrapped samples.
%   .wReplace = boolean, Default=1=bootstrap sample with-replacement. 0=sample
%               without-replacement.
%   .testFrac = fraction of data to be use for early stopping or indirect 
%               fitting.
%   .randSeed = random seed used to initialize the state of rand.
%   .optimOpt = options required by the specified optimization algorithm.
%               Default uses trnGradDesc.m for training.
% OUTPUT:
% [strfArray] = structure array of structure containing [options.nResamp] 
%               fitted models.
%  [optArray] = structure array of options
%   [bootIdx] = Index of which samples were used in each bootstrap sampling.
%
% SEE ALSO: resampJackknife, resampCrossVal
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jul 2007)
%
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='resampBootstrap';
optRange.funcName={'resampBootstrap'};

optDeflt.nResamp=10;
optRange.nResamp=[1,10000];
optDeflt.wReplace=1;
optRange.wReplace=[0,1];
optDeflt.testFrac=0.10;
optRange.testFrac=[0,1];
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
    warning(['You are using a fitting algorithm that uses an early stopping set, and have testFrac=0, therefore the testFrac has been set to 0.10']);
    options.testFrac=0.10;
    pause(2);
elseif options.testFrac>0 & (~isfield(options.optimOpt, 'earlyStop') | options.optimOpt.earlyStop == 0);
    warning(['You are using a fitting algorithm that does not use an early stopping set, and have testFrac>0, therefore some' ...
        ' data is not being ignored from the training set']);
    pause(2);
end

% Get sample size & Check Group Index consistency
% --------------------
groupIdx=globDat.groupIdx;
sampSiz=length(datIdx);
if ~isempty(groupIdx)
  groupLen=length(groupIdx);
  if groupLen~=sampSiz
    error('resampBootstrap >< group index must have same length as data index');
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
gTestLen=round(options.testFrac*nGroup);
testLen=round(options.testFrac*sampSiz);
trainLen=sampSiz-testLen;


% Compute models over bootstrapped samples
% --------------------
rand('state',options.randSeed);
for ii=1:options.nResamp
  disp(sprintf('----- Bootstrap Sample %d/%d',ii,options.nResamp));
  rsIdx=randperm(nGroup);
  tstIdxSet=findIdx(rsIdx(1:gTestLen),groupIdx);
  trnIdxSet=findIdx(rsIdx((gTestLen+1):end),groupIdx);
  
  if options.wReplace
    trainIdx=randsample(trnIdxSet,trainLen,true);
    testIdx=randsample(tstIdxSet,testLen,true);
  else
    trainIdx=trnIdxSet;
    testIdx=tstIdxSet;
  end
  
  [strfArray(ii),optArray(ii)]=feval(options.optimOpt.funcName,strf, ...
    datIdx(trainIdx),options.optimOpt,datIdx(testIdx));
  
  bootIdx(:,ii)=testIdx;
end


