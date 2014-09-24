%% By making some reasonable assumptions and using iterative gradient
%% descent methods, we have pushed RTA and RTC about as far as possible.
%% Other, more complex models may be desireable, however, especially when
%% attempting to model areas beyond V1. Currently our best model consists
%% of preprocessing the stimuli with a pyramid of 3D Gabor filters, that
%% vary in scale, spatial frequency, temporal frequency, orientation and 
%% position which are then either rectified or combined and squared, to give
%% an assortment non-linear filters resembling of simple and complex cells.
%% These filters will function as the linearizing transform applied to the
%% stimulus. The fitting routine determines the linear combination of these 
%% filters that best fits the data. There are many free parameters in this 
%% model that determine all aspects of the set of Gabor filters and the 
%% non-linearity. The documentation in preprocWavelets3d explains what 
%% each of these free parameters do.


%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_5_linearization');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))

global globDat


%% load a 20x20x20000 natural movie
movPath = fullfile(strflabDir, 'fakedata/mov.mat');
load(movPath);

%% let's take the whole movie and later separate the data into training and
%% validation sets, so we can see how well our models predict
rawStim = single(mov);

%% Lets make a complex cell with a preferred direction and a certain temporal
%% frequency
gparams = [0.5000 0.5000 0 5.8654 1.4205 0.0852 0.2816 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
%% You should probably look at each of these filters to see what results to
%% expect, i.e. imagesc(gabor(:,:,3))

gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
rawStim = reshape(rawStim, [10*10 20000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);
resp = [zeros([4 1]); resp(1:end-4)];

%% Now let's add noise to the computed response. 
resp = resp + 1.5*std(resp)*randn(size(resp));

%% delay the response as before
resp = [zeros([4 1]); resp(1:end-4)];

%% Since it is possible to have more than one spike per frame, let's try
%% a little more information preserving threshhold, with a high threshhold 
%% for 2 spikes and a little bit lower threshhold for 1 spike
resp(resp<75) = 0;
resp(resp>=150) = 2;
resp(resp>=75) = 1;

%% Set the params to the defaults for preprocWavelets3d. Modify and 
%% experiment to see how the size of stim changes. Each column of stim
%% represents the response of a filter to the stimulus
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
options.stepSize = 1e-04;

%% We'll try using 80% of the data for the training set
trainingIdx = [1:floor(.8*globDat.nSample)];

%% And the remaining 20% for the stopping set
stoppingIdx = [floor(.8*globDat.nSample)+1:globDat.nSample];

%% Train and viusalize the STRF
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);

%% Now let's globalize our validation data and use the trained
%% STRF to predict the response
strfData(stimVal,respVal)
testingIdx = [1:globDat.nSample];
[strfTrained,predresp]=strfFwd(strfTrained,testingIdx);

%% Let's remove any NaN's from the prediction and look at the
%% correlation between the prediciton and actual response
nonnanidx = find(~isnan(predresp));
predcorr1 = corr(predresp(nonnanidx),respVal(nonnanidx))

%% That's pretty good for as noisy as we made the response.
%% Now let's graph the predicion and response to compare:
figure; plot(respVal(nonnanidx), 'r'); hold on; plot(predresp(nonnanidx), 'b')

%% Clearly the prediction and actual response are of different kinds.
%% The prediction of the linear model can be thought of as an average
%% spike rate, i.e. what would result from averaging several trials of
%% spike responses together. It will of course look different than
%% the spikes on a single trial, but we seem to have captured the major
%% increases and decreases in spike rate through the trial quite well.

%% Now lets see what we have as far as filter estimates
preprocWavelets3dVis(strfTrained);

%% We see that the final filter qualitatively looks decent and captures
%% filter characteristcs well, but is pretty noisy because the final
%% model includes a weight on every stimulus channel. This is because 
%% gradient descent computes the gradent for all the channels and takes a
%% small step down the entire gradient, which usually results in all terms
%% having some weight. If however we only take a step down the steepest 
%% direction of the gradient, and use early stopping, we should end up
%% with a much sparser answer. The practice of only taking steps along 
%% the steepest direction of the gradient at each iteration is known as 
%% coordinate descent. As discussed in the STRFlab paper, using coordinate 
%% descent and early stopping cooresponds to a sparse prior over the 
%% parameter values. Because each of our stimulus channels corresponds
%% to an entire non-linear filter, rather than a single pixel or covariance
%% of 2 pixels, it is very reasonable to want a STRF that describes the 
%% system using as few terms as possible.

%% Put the training data back as the global variable
strfData(stim,resp)

%% Coordinate descent is an option in the trnGradDesc routine and can be 
%% enabled:
options.coorDesc = 1;

%% I've chosen a larger step size to speed things up. Again this may need
%% to be adjusted for your particular situation
options.stepSize = 1e-03;

%% Now let's train and viusalize the STRF fit with coordinate descent
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);

%% Now let's globalize our validation data and use the trained
%% STRF to predict the response
strfData(stimVal,respVal)
testingIdx = [1:globDat.nSample];
[strfTrained,predresp]=strfFwd(strfTrained,testingIdx);

%% Let's remove any NaN's from the prediction and look at the
%% correlation between the prediciton and actual response
nonnanidx = find(~isnan(predresp));
predcorr2 = corr(predresp(nonnanidx),respVal(nonnanidx))

%% That's better than before, and we used fewer parameters to get a better
%% prediction. Clearly coordinate descent is the preferred method here.

%% Now let's graph the predicion and response to compare:
figure; plot(respVal(nonnanidx), 'r'); hold on; plot(predresp(nonnanidx), 'b')

%% Again we seem to have captured the major increases and decreases in spike 
%% rate through the trial quite well.

%% Now lets see what we have as far as filter estimates
preprocWavelets3dVis(strfTrained);

%% This is obviously a much more sparse answer than we obtained previously.
%% The sparseness makes it more clear what is going on, and this strf gives
%% better predictions, so it is obviously preferable.