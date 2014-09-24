
%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_6_resampling');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))

global globDat

%% load a 20x20x20000 natural movie
movPath = fullfile(strflabDir, 'fakedata/mov.mat');
load(movPath);

%% let's only take a part of it so it's directly comparable to above
rawStim = single(mov);

%% Let's create a complex gabor filter and the response as before
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
rawStim = reshape(rawStim, [10*10 20000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);

%% add some noise to the computed response.
resp = resp + 1.5*std(resp)*randn(size(resp));

%% delay the response as before
resp = [zeros([4 1]); resp(1:end-4)];

%% Now we're going to malicioiusly corrupt some of the data. In this worst
%% case scenario the corrupted data does not look obviously different than
%% the real data and would be difficult to impossible to locate and remove
resp(10001:15000) = mean(resp) + std(resp)*(randn(5000,1));

%% Set the spiking thresholds as before
resp(resp<150) = 0;
resp(resp>=300) = 2;
resp(resp>=150) = 1;

%% Set the params to the defaults for preprocWavelets3d. 
rawStim = reshape(rawStim, [10 10 20000]);
params = preprocWavelets3d;
[stim,params] = preprocWavelets3d(rawStim, params);

%% Now let's break the data into Validation and training sets
%% Validation:
stimVal = stim(15001:end,:);
respVal = resp(15001:end);

%% Training:
stim = stim(1:15000,:);
resp = resp(1:15000);

%% Globalize the stim and resp
strfData(stim,resp)

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% set the options to defaults for trnGradDesc
options = trnGradDesc;

%% turn on early stopping
options.earlyStop = 1;
options.display = -1;
options.nDispTruncate = 0;
options.coorDesc = 1;

%% I've chosen a larger step size to speed things up. Again this may need
%% to be adjusted for your particular situation
options.stepSize = 1e-03;

%% We'll try using 80% of the data for the training set
trainingIdx = [1:floor(.8*globDat.nSample)];

%% And the remaining 20% for the stopping set
stoppingIdx = [floor(.8*globDat.nSample)+1:globDat.nSample];

%% Train and viusalize the STRF
strfTrained1=strfOpt(strf,trainingIdx,options,stoppingIdx);

%% Now let's globalize our validation data and use the trained
%% STRF to predict the response
strfData(stimVal,respVal)
testingIdx = [1:globDat.nSample];
[strfTrained1,predresp]=strfFwd(strfTrained1,testingIdx);

%% Let's remove any NaN's from the prediction and look at the
%% correlation between the prediciton and actual response
nonnanidx = find(~isnan(predresp));
predcorr1 = corr(predresp(nonnanidx),respVal(nonnanidx))



options = resampBootstrap;
options.optimOpt = trnGradDesc;
options.optimOpt.coorDesc = 1;
options.optimOpt.earlyStop = 1;
options.optimOpt.nDispTruncate = 0;
options.optimOpt.display = -1;
options.optimOpt.stepSize =  1e-03;
options.testFrac = 0.20;
options.nResamp = 5;

strfData(stim,resp)

trainingIdx = [1:globDat.nSample];
strfTrained_tmp=strfOpt(strf,trainingIdx,options);
strfTrained2 = strfTrained_tmp(1);
strfTrained2.b1 = mean(cat(2, strfTrained_tmp.b1), 2);
strfTrained2.w1 = mean(cat(4, strfTrained_tmp.w1), 4);

strfData(stimVal,respVal)
testingIdx = [1:globDat.nSample];
[strfTrained2,predresp]=strfFwd(strfTrained2,testingIdx);

%% Let's remove any NaN's from the prediction and look at the
%% correlation between the prediciton and actual response
nonnanidx = find(~isnan(predresp));
predcorr2 = corr(predresp(nonnanidx),respVal(nonnanidx))