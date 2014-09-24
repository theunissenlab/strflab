%Typical vision setup (Not part of strflab)
%---------------------------------------
% load /auto/k5/prenger/strflab/fakedata/sampledata_nonan.mat
% rmpath(genpath('/auto/k2/share/strflab'));  % remove current STRF lab path
% addpath(genpath('/auto/k2/share/strflabGO'));  % add new STRF lab path

% stimEst=x(1:900,:);
% respEst=y(1:900);
% stimVal=x(901:1000,:);
% respVal=y(901:1000);
% load /auto/k1/moliver/code/midtestdata.mat

x=randn(10000,100);

h0=randn(100,1);
y=x*h0;


z=find(y<=0);
y(z)=0;
z=find(y>0);
y(z)=1;


stimEst=x(1:9000,:);
respEst=y(1:9000);
stimVal=x(9001:10000,:);
respVal=y(9001:10000);


% Declaring global variables.
%-----------------------------------
global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(stimEst,respEst);


%Initialize strf
%--------------------------------------
% strf=midInit(100,[0:1]);
strf=midInit(100,[0:1],'nonparametric',5,0.2,'const');


%Train with Gradient Descent with early stopping
%----------------------------------
% options=trnGradDesc
% options.display=-5;
% options.earlyStop=1;
% options.stepSize = .005;
% options.adaptive = 0;
% options.coorDesc = 1;
% trainingIdx=1:floor(9*globDat.nSample/10);  % generate index of training samples.
% earlyStopIdx=(trainingIdx(end)+1):globDat.nSample;  % generate index of early stopping samples.
% strfTrained=strfOpt(strf,trainingIdx,options,earlyStopIdx);

% 
% %Train using Jackknife resampling
% %----------------------------------
% options=resampJackknife;
% options.nResamp = 3;
% options.optimOpt=trnGradDesc;
% options.optimOpt.coorDesc=1;
% options.optimOpt.display=-5;
% options.optimOpt.earlyStop=1;
% options.optimOpt.stepSize = .005;
% options.optimOpt.adaptive = 0;
% 
% trainingIdx=1:globDat.nSample;  % generate index of estimation samples.
% strfArray=strfOpt(strf,trainingIdx,options);
% strfTrained = strfArray(1);
% strfTrained.w1 = mean(cat(4, strfArray.w1), 4);
% strfTrained.b1 = mean(cat(2, strfArray.b1), 2);
% 
% keyboard
% Train with plain SCG
% ----------------------------------
options=trnSCG;
options.display=-5;

trainingIdx=1:globDat.nSample;  % generate index of estimation samples.
strfTrained=strfOpt(strf,trainingIdx,options);


% Prediction: **** LOOK AT THIS ****
%----------------------------------
strfData(stimVal,respVal);  % Set the validation data to be the global variable, so you don't predict on estimation data.
[strfTrained,predResp]=strfFwd(strfTrained,1:globDat.nSample);
figure; plot(predResp,respVal);

nonNanIdx=intersect(find(~isnan(predResp)),find(~isnan(respVal)));
corr2(predResp(nonNanIdx),respVal(nonNanIdx))


%rmpath(genpath('/auto/k2/share/strflabind'));  % remove new STRF lab path
%addpath(genpath('/auto/k2/share/strflab'));  % add back current STRF lab path
