%Typical auditory setup (Not part of strflab)
%---------------------------------------
load /auto/k2/share/strflabGOLD/fakedata/sampleDataAuditory.mat
addpath(genpath('/auto/k2/share/strflabGOLD'));  % add new STRF lab path
% 
% stimEst=x(1:900,:);
% respEst=y(1:900);
% stimVal=x(901:1000,:);
% respVal=y(901:1000);
estCutoff = 35000;
 x = x - mean(x,2)*ones(1,size(x,2));
 x = x ./ (std(x')'*ones(1,size(x,2)));
stimEst = (x(:,1:estCutoff)');
respEst = mean(t(:,1:estCutoff))';
stimVal = (x(:,(1+estCutoff):end)');
respVal= mean(t(:,(1+estCutoff):end))';
% Declaring global variables.
%-----------------------------------
global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(stimEst,respEst);


%Initialize strf
%--------------------------------------
if ~exist('strf','var')
strf=linInit(size(stimEst,2),[0:40]);
end


%Train with Gradient Descent with early stopping
%----------------------------------

options=trnGradDesc;
options.display=-5;
options.earlyStop=1;
options.adaptive = 0;
options.stepSize = 2e-7;%.05; % .000001;%  
options.funcName ='trnGradDesc';
options.maxIter = 1000;
trainingIdx=1:floor(9*globDat.nSample/10);  % generate index of training samples.
earlyStopIdx=(trainingIdx(end)+1):globDat.nSample;  % generate index of early stopping samples.
globDat.stim = globDat.stim - mean(mean(globDat.stim));
globDat.resp = globDat.resp - mean(mean(globDat.resp));
tic;
[strfTrained,options]=strfOpt(strf,trainingIdx,options,earlyStopIdx);
toc;
foundStepSize = options.stepSize;
if all(strfTrained.w1(:) == strf.w1(:))
    disp(['The first try didn''t optimize the STRF at all.  Maybe the step size was too big.  Trying again...']);
tic;
[strfTrained,options]=strfOpt(strfTrained,trainingIdx,options,earlyStopIdx);
toc;
end

%{
%Train using Jackknife resampling
%----------------------------------
options=resampJackknife;
options.nResamp = 3;
options.optimOpt.coorDesc=1;
options.optimOpt.display=-5;
options.optimOpt.earlyStop=1;
if exist('foundStepSize','var')
    options.optimOpt.stepSize = foundStepSize;
else
    options.optimOpt.stepSize = .05;
end
options.optimOpt.funcName ='trnGradDescAdaptive';
trainingIdx=1:globDat.nSample;  % generate index of estimation samples.
tic;
[strfTrained,options]=strfOpt(strf,trainingIdx,options);
toc;


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
plot(predResp,respVal);

nonNanIdx=intersect(find(~isnan(predResp)),find(~isnan(respVal)));
corr2(predResp(nonNanIdx),respVal(nonNanIdx))


%rmpath(genpath('/auto/k2/share/strflabind'));  % remove new STRF lab path
%addpath(genpath('/auto/k2/share/strflab'));  % add back current STRF lab path
%}