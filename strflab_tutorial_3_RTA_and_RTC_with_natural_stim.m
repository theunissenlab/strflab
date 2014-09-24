%% Now let's try to fit our RTA and RTC models using natural movies

%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_3_RTA_and_RTC_with_natural_stim');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))

global globDat

%% load a 20x20x20000 natural movie
movPath = fullfile(strflabDir, 'fakedata/mov.mat');
load(movPath);

%% let's only take a part of it so it's directly comparable to above
rawStim = single(mov(1:10,1:10,1:15000));

%% let's create some fake data for a simple cell with a 2D Gabor filter
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 1], gparams);
gabor = reshape(gabor, [10*10 1]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp = [zeros([4 1]); resp(1:end-4)];

%% reshape and preprocess the stimuli
rawStim = reshape(rawStim, [10 10 15000]);
params.RTAC = [1 0];
[stim,params] = preprocRTAC(rawStim, params);

%% create a new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% globalize the data
strfData(stim,resp)

%% set options and index as before
%% notice we increased the maximum number of iterations, due to correlations
%% in the natural stimulus, the training will take longer to converge 
options=trnSCG;
options.display=-5;
options.maxIter = 300;
trainingIdx = [1:globDat.nSample];

%% Train and visualize the strf
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% It should be obvious that the filter has again been very well recovered by
%% RTA. Doing the RTA in terms of a linear regression is actually equivalent
%% to doing a whitened RTA, which involves inverting the stimulus covariance
%% matrix to remove stimulus correlations. People generally used whitened
%% stimuli to avoid having to invert the covariance matrix. Doing the linear 
%% regression by gradient descent (and quickly descending all the way using SCG)
%% avoids having to invert the covariance matrix.
%% The advantages of STRFlab should start becoming aparent to you...

%% Now lets make some more realistic fake data as before, using a 3D Gabor filter
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp = [zeros([4 1]); resp(1:end-4)];

%% If you plot the response you might be surprised at how much larger the 
%% magnitude is with natural images compared to the Gaussian white noise
%% used earlier. Gabor filters, however, are a great basis set for representing
%% natural images so it shouldn't be so surprising that natural movies have
%% large projections onto gabors. 
%% But for our fake data, we will need to set a much higher threshhold to get 
%% a reasonable density of spikes
resp(resp<50) = 0;
resp(resp>=50) = 1;

%% reshape and preprocess the stimuli
rawStim = reshape(rawStim, [10 10 15000]);
params.RTAC = [1 0];
[stim,params] = preprocRTAC(rawStim, params);

%% Put the new response into the global variable
strfData(stim,resp);

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% Now train and visualize the new STRF using the same options as before
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% Again we see the filter has been reasonably well recovered. It is noiser
%% than before due to the information lost using the hard spiking threshhold
%% and limited data.

%% Now lets try a complex cell with natural movies
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);
resp = [zeros([4 1]); resp(1:end-4)];

%% Since it is possible to have more than one spike per frame, let's try
%% a little more information preserving threshhold, with a high threshhold 
%% for 2 spikes and a little bit lower threshhold for 1 spike
resp(resp<150) = 0;
resp(resp>=300) = 2;
resp(resp>=150) = 1;

%% To model a complex cell we turn on RTC in addition to RTA
params = preprocRTAC;
params.RTAC = [1 1];
params.locality = 1;

%% Make sure rawStim is the proper size before sending to preprocessing
rawStim = reshape(rawStim, [10 10 15000]);
[stim,params] = preprocRTAC(rawStim, params);
strfData(stim,resp)

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% Train and visualize the model
options=trnSCG;
options.display=-5;
options.maxIter = 1500;
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% We see that the 1st principle components of the RTC matrix around 
%% the peak delay recover the filter as well as in the white noise case. 
%% By framing the RTC problem in terms of regression and using gradient
%% descent, we were able to perform a RTC analysis with natural stimuli 
%% without inverting a matrix of covariance-of-the-covariance, 4th order 
%% terms. This is clearly a powerful method for performing a whitened STC

%% Note that the RTA estimate appears to have some structure. This is due
%% to correlations in natural stimuli beyond what can be accounted for by
%% inverting the stimulus covariance matrix. The filter we created was 
%% essentially a polarity insensitive edge detector. Therefore averaging 
%% over the natural stimuli that modulate the filter often results in a
%% bright blob centerd on the receptive field due to biases in the stimulus. 