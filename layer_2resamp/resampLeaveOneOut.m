function [strfArray,optArray,cvResult]=resampLeaveOneOutCV(strf,datIdx,options)
%function [strfArray,optArray,cvResult]=resampCrossVal(strf,datIdx,options)
%
% Resample the data with replacement and compute model STRF parameters for
% each set of resampled data.
%
% INPUT:
%      [strf] = model structure obtained via upper level *Init functions
%    [datIdx] = a vector containing indices of the samples to be use in the 
%               resampling fitting.
%   [options] = option structure containing fields:
%   .funcName = char array, name of this function.
%      .nFold = # of even splits of the data into smaller subsets. Each time
%               1 of the subset is used as test data and trained on the rest.
%               The prediction is generated from these test data.
%    .nResamp = # of cross validation rounds. Each round fits 1 model and
%               predict on the hold-out test data.
%  .randomize = boolean, Default=1=randomize the order of the samples first
%               before splitting it n-fold. 0=Keep the original data order.
%   .randSeed = random seed used to initialize the state of rand.
%   .optimOpt = options required by the specified optimization algorithm.
%               Default uses trnGradDesc.m for training.
% OUTPUT:
% [strfArray] = structure array of structure containing [options.nFold] 
%               fitted models.
%  [optArray] = structure array of options
%  [cvResult] = structure array of CV results containing fields:
%       .pred = model predicted response from each set of test data.
%       .true = true measure response in the test data.
%    .testIdx = index of test data samples.
%
% SEE ALSO: resampJackknife, resampBootstrap
%
% By Michael Wu  --  waftingpetal@yahoo.com (Jul 2007)
% Bug Fix to allow more than 5-fold resampling by Michael Oliver (Jul 2009)
% Bug Fix to properly return each CV prediction by Michael Oliver (Jul 2009)
%
% ====================


% Set default option values
% --------------------
global globDat;
optDeflt.funcName='resampLeaveOneOut';
optRange.funcName={'resampLeaveOneOut'};
optDeflt.testFrac=0.1;
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


% Get sample size & Check Group Index consistency
% --------------------
groupIdx=globDat.groupIdx;
sampSiz=length(datIdx);
if ~isempty(groupIdx)
  groupLen=length(groupIdx);
  if groupLen~=sampSiz
    error('resampCrossVal >< group index must have same length as data index');
  end
  uniqGS=unique(groupIdx);
  nGroup=length(uniqGS);
else
  uniqGS=1:sampSiz;
  groupIdx=uniqGS;
  nGroup=sampSiz;
end

  rsIdx=uniqGS;



% Compute models & prediction over nFold split of data
% --------------------
for ii=1:nGroup
  disp(sprintf('----- %d Fold Cross Validation: %d/%d ',nGroup,ii,nGroup));
  gTrainIdx=setdiff(rsIdx,ii);
  testIdx=findIdx(ii,groupIdx);
  trainIdx=findIdx(gTrainIdx,groupIdx);
  
  if options.testFrac > 0
    esIdx=randsample(trainIdx,floor(options.testFrac*size(trainIdx,2)),false);
    trainIdx = setdiff(trainIdx, esIdx);
  else
    esIdx=testIdx;
  end
  
  [strfArray(ii),optArray(ii)]=feval(options.optimOpt.funcName,strf, ...
    datIdx(trainIdx),options.optimOpt,datIdx(esIdx));

  % Store results: true & predicted response, and samples index of test data
  [strfArray(ii), cvResult(ii).pred]=strfFwd(strfArray(ii),datIdx(testIdx));
  cvResult(ii).true=globDat.resp(testIdx,:);
  cvResult(ii).testIdx=testIdx;
end


