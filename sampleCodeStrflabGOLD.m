%Typical vision setup (Not part of strflab)
%---------------------------------------
addpath(genpath('/auto/k2/share/strflabGOLD'));  % add new STRF lab path
load /auto/k2/share/strflabGOLD/fakedata/sampledata_nonan.mat

stimEst=x(1:900,:);
respEst=y(1:900);
stimVal=x(901:1000,:);
respVal=y(901:1000);


% Declaring global variables.
%-----------------------------------
global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(stimEst,respEst);


%Initialize strf
%--------------------------------------
strf=linInit(100,[0:4]);


%Train with Gradient Descent with early stopping
%----------------------------------
options=trnGradDesc
options.display=-5;
options.earlyStop=1;
options.stepSize = .05;
options.adaptive = 0;

trainingIdx=1:floor(9*globDat.nSample/10);  % generate index of training samples.
earlyStopIdx=(trainingIdx(end)+1):globDat.nSample;  % generate index of early stopping samples.
strfTrained=strfOpt(strf,trainingIdx,options,earlyStopIdx);


%Train using Jackknife resampling
%----------------------------------
options=resampJackknife;
options.nResamp = 3;
options.optimOpt=trnGradDesc;
options.optimOpt.coorDesc=1;
options.optimOpt.display=-5;
options.optimOpt.earlyStop=1;
options.optimOpt.stepSize = .05;
options.optimOpt.adaptive = 0;

trainingIdx=1:globDat.nSample;  % generate index of estimation samples.
strfTrained=strfOpt(strf,trainingIdx,options);


%Train with plain SCG
%----------------------------------
options=trnSCG;
options.display=-5;

trainingIdx=1:globDat.nSample;  % generate index of estimation samples.
strfTrained=strfOpt(strf,trainingIdx,options);


% Prediction: **** LOOK AT THIS ****
%----------------------------------
strfData(stimVal,respVal);  % Set the validation data to be the global variable, so you don't predict on estimation data.
[strfTrained,predResp]=strfFwd(strfTrained,1:globDat.nSample);
plot(predResp,respVal, '.');

nonNanIdx=intersect(find(~isnan(predResp)),find(~isnan(respVal)));
corr2(predResp(nonNanIdx),respVal(nonNanIdx))


%rmpath(genpath('/auto/k2/share/strflabind'));  % remove new STRF lab path
%addpath(genpath('/auto/k2/share/strflab'));  % add back current STRF lab path
