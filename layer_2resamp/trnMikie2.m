function [strf]=trnMikie2(strf,stim,resp,options)
%function [strf]=trnMikie2(strf,stim,resp,options)
%
% Sample training algorithm for Mikie. This training algorithm will
% use gradient descent as the low level optimizer to train a model.
% It will Jackknife the data and re-estimate the models over the 
% Jackknife samples. And model selection is done by model averaging.
%




% Set default option values
% --------------------
optDeflt=resampJackknife;
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


% Jackknife the data and estimate models
% --------------------
[strfArray]=resampJackknife(strf,stim,resp,options);
parmArray=strfPak(strfArray);


% Model selection by averaging the model parameters.
% --------------------
avgMod=mean(parmArray,1);
strf=strfUnpak(strfArray,avgMod);


