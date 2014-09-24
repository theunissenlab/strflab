function [strf]=trnMikie(strf,stim,resp,options)
%function [strf]=trnMikie(strf,stim,resp,options)
%
% Sample training algorithm for Mikie. This training algorithm will
% use gradient descent as the low level optimizer to train a model.
% The model have 1 hyper-parameter. And its value is selected by
% 3-fold cross validation. And the model selection metric is 
% the correlation coefficient. 




% Set default option values
% --------------------
optDeflt=resampCrossVal;
optDeflt.optimAlg='trnGradDesc';
optRange.optimAlg={'trnGradDesc','trnSCG','trnDirectFit', ...
  'trnBoosting','trnSimAnneal','trnQuadProg'};

if nargin<4
  options=optDeflt;
else
  options=defaultOpt(options,optDeflt,optRange);
end

if nargin<1
  strf=optDeflt;
  return;
end


% Choose a list of hyper-parameters.
% --------------------
hyperParm=logspace(log(1),log(1000),10);
nHyperParm=length(hyperParm);


% Loop over hyperparameters.
% --------------------
for ii=1:nHyperParm
  strf.hyperParam=hyperParm(ii);  % set the hyper-param in the model structure.
  [strfArray,cvResult]=resampCrossVal(strf,stim,resp,options);  % do CV
  
  % use correlation coefficient as the validation metric
  for jj=1:options.nFold
    predScore(jj)=corr2(cvResult(jj).pred,cvResult(jj).true);
  end
  r(ii)=mean(predScore);  % get a mean corr coeff

end


% Model selection by picking the highest corr coeff.
% --------------------
[rMax,maxIdx]=max(r);

% Set the hyper-param again as the best hyper-param from CV and train the model again
% --------------------
strf.hyperParam=hyperParm(maxIdx);
strf=feval(options.optimAlg,strf,stim,resp,option.trnOpt);

