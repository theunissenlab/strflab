%%%  STRFlab Tutorial - version 1.0  %%%

%% This tutorial is designed to familiarize you with using all of the 
%% features of STRFlab. STRFlab need only be in your path and you should
%% be able to copy and paste all the code here and have it run.

%% We begin with some simple examples of the response triggered average (RTA) 
%% and response triggered covariance (RTC) reformulated in terms of a linear
%% regression, which can be solved using iterative error minimization 
%% techniques and thus take advantage of STRFlab's architecture.
%% In this formulation RTA is equivalent to finding a weight for each pixel
%% (via gradient descent) such that the weighted sum gives the predicted
%% response; and RTC is equivalent to finding weights for each pair-wise
%% product of pixels (the stimulus covariance matrix) such that the weighted 
%% sum gives the predicted response. PCA must be performed on the RTC
%% weights to recover the filters (the visualization routine does this
%% automatically)

%% Make sure you have STRFlab in your path
addpath(genpath('/auto/k2/share/strflabGOLD'))

%% Let's start by making some fake data

%% First we make Gaussian random noise stimuli
rawStim = randn([10 10 15000]);

%% Now we define a filter using make3dgabor which yields an even and odd
%% phase vertical Gabor filter
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 1], gparams);

%% We are going to start using only the even filter to model a simple cell
%% So, first reshape the 2-D filter and stimulus into sets of 1-D vectors
gabor = reshape(gabor, [10*10 1]);
rawStim = reshape(rawStim, [10*10 15000]);

%% dotdelay computes the filter response of the stimulus
resp = dotdelay(gabor, rawStim);

%% To make the example more realistic we insert a delay between the
%% stimulus and response
resp = [zeros([4 1]); resp(1:end-4)];

%% We reshape the stimulus before preprocessing because all stimuli 
%% preprocessing routines are designed to work with the raw stimuli as
%% 3D (or 4D for color) matrices
rawStim = reshape(rawStim, [10 10 15000]);

%% Now we are done making fake data. We are at the point of just having
%% the raw stimulus and a set of responses as in a real experiment and 
%% can now proceed normally

%% There is a stimulus preprocessing routine for doing RTA and RTC,
%% called preprocRTAC. The parameter params.RTAC determines whether to
%% compute RTA or RTC. To do only RTA set it to [1 0], to do RTC only set 
%% to [0 1], or to do both set to [1 1]. We start with just RTA. In the 
%% case of RTA, the preprocessing routine doesn't have to do much since
%% there is no transformation of the stimulus, it just formats the
%% stimulus to work with STRFlab
params.RTAC = [1 0];
[stim,params] = preprocRTAC(rawStim, params);

%% Now we create a STRF structure variable that contains all the model
%% parameters that will be fit to the data. We create a Generalized 
%% Linear Model (GLM) with the linInit function. The first argument is 
%% the number of parameters in the fitting, and the second argument is 
%% the number of delays to consider between stimulus and response
strf = linInit(size(stim,2), [0:8]);

%% Set the bias term of the GLM to the mean response. This generally
%% speeds up the fitting process
strf.b1 = mean(resp);

%% Put the preprocessing parameters into the strf structure to keep a 
%% record of what the fitted model represents and to allow the 
%% visualization routines to work properly
strf.params = params;

%% STRFlab uses a global variable to share data between functions.
%% This must be instantiated here
global globDat

%% strfData puts the stimulus and response that will be used to fit the
%% model into the globab variable
strfData(stim,resp)

%% To train this model we are going to use the scaled conjugate gradient
%% technique in the function trnSCG. Calling it with no arguments returns
%% the default set of options, which can then be edited. We will change it
%% to graphically display every five iterations and only perform at most
%% 100 iterations
options=trnSCG;
options.display=-5;
options.maxIter = 100;

%% Create an index of the samples used to train the model. 
trainingIdx = [1:globDat.nSample];

%% strfOpt is a general routine that takes any strf structure and set of
%% options and uses the appropriate error and fitting functions. It returns
%% a new structure with the model parameters fit to the data. Because of
%% the display option above, this will open a window showing model 
%% parameter values and error values
strfTrained=strfOpt(strf,trainingIdx,options);

%% Now we can use the visualization routine that corresponds to the preproc
%% routine to visualize the STRF
preprocRTAC_vis(strfTrained);

%% The visualization should reveal that RTA has recovered the filter at 
%% the 4th delay as intended.

%% Now lets consider a little more complicated filter; one that has a time-
%% course of 5 frames and then make the fake response data
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp = [zeros([4 1]); resp(1:end-4)];

%% Let's further make this more realistic by postively rectifying the filter
%% response and setting a spiking threshhold. Now this response data closely
%% resembles data from a neuron. The threshhold set here for this fake data 
%% is arbitrary and chosen to give a resonable density of spikes
resp(resp<3) = 0;
resp(resp>=3) = 1;

%% Put the new response into the global variable
strfData(stim,resp);

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% Now train and visualize the new STRF
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% We see that even with a simulated spiking neuron we can recover the filter
%% with RTA

%% The preceding showed how the filter for a simple cell could be recovered
%% using RTA. Let's now consider a complex cell. We use both the even and odd
%% Gabor filters obtained previously. The sum of the squared responses of each
%% filter corresponds to the energy model of a complex cell. We then create a
%% delay, rectify and threshhold as before.
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);
resp = [zeros([4 1]); resp(1:end-4)];
resp(resp<4) = 0;
resp(resp>=4) = 1;

%% To model a complex cell we turn on RTC in addition to RTA
params.RTAC = [1 1];

%% The parameter covtime specifies whether to compute the covariance of the
%% stimulus over time. Setting it to 0 will calculate the covariace across 
%% each frame separately, without calculating covariace across time
params.covtime = 0;

%% In order to constrain the model to a reasonable number of parameters
%% (unlike standard RTC) we can consider only local pairwise statistics.
%% The locality parameter determines how many pixels away from each pixel to
%% include pairwise statistics for. Setting to to one includes all bordering
%% pixels.
params.locality = 1;

%% Make sure rawStim is the proper size before sending to preprocessing
rawStim = reshape(rawStim, [10 10 15000]);

%% The preprocessing routine is now going to return the raw stimulus, to
%% compute the RTA, and a limited portion (determined by the locality 
%% paramter) of the covariance matrix of the stimulus, to compute the RTC,
%% all formatted for STRFlab
[stim,params] = preprocRTAC(rawStim, params);
strfData(stim,resp)

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% Train and visualize the model
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% It should be obvious that the filter estimate obtained by RTC for this model 
%% complex cell is far noisier than those for model simple cells obtained with
%% RTA. This is largely because we are trying to fit a large model with limited
%% data. I leave it to the user to experiment with increasing the size of the 
%% stimulus set (i.e. to 30000 frames), fitting only the RTC (params.RTAC = [0 1]),
%% including more pairwise statistics (i.e. params.locality = 2), fitting the
%% strf only at the peak delay (i.e. strf = linInit(size(stim,2), [4]) )
%% This problem of model size for RTC will only get worse as the number of pixels
%% in the stimulus increases; so far we have only considered a very small 10X10
%% stimulus, but for real experiments one would probably prefer something larger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now let's try to fit our RTA and RTC models using natural movies

%% load a 20x20x20000 natural movie
load /auto/k2/share/strflabGOLD/fakedata/mov.mat

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
options.maxIter = 250;
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
params.RTAC = [1 1];
params.locality = 1;
params.covtime = 0;

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

%% Yuck! Isn't that ugly. We see something like our desired filter at 
%% the peak in the first set of RTC filters, but its pretty gross.
%% Again we have many parameters, limited data made noisy by
%% threshholding, and lots of stimulus correlations in natural movies. 
%% But lets try something...
options.maxIter = 15;
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% Lo and behold, the filters look much cleaner. What is going on here?
%% We didn't even descend the error gradient all the way, and yet it 
%% looks better?!?! This is called early-stopping and can prevent over-
%% fitting to noise in the data and stimulus correlations. Gradient 
%% descent fits the parameters that account for the most variance first,
%% and since most of the variance in natural images is low frequency, 
%% early stopping is effectively low pass filtering the result, making
%% it look better. Can we do early stopping in a principled way rather 
%% than just picking an arbitrary number of iterations? Of course, or I 
%% wouldn't have brought it up. trnSCG tries to descend the gradient as 
%% quickly as possibletaking steps of various sizes and is not well suited 
%% to doing early stopping. So we turn now to another fitting routine 
%% trnGradDesc. This function can take an index of data samples to use as 
%% a stopping set. The gradient will be computed on the training set as 
%% usual and a small step of a specified size will be taken down the 
%% gradient. The error is also assessed at each step on the stopping set 
%% and if the error goes up too much on this separate set, the training 
%% will stop.

%% set the options to defaults for trnGradDesc
options = trnGradDesc;

%% turn on early stopping
options.earlyStop = 1;

%% graphically display every iteration
options.display = -1;

%% don't truncate first steps
options.nDispTruncate = 0;

%% In practice the step size for trnGradDesc may need to be tuned for 
%% each case. To do so just set a step size to something small enough 
%% error doesn't occilate on the training set with each step and doesn't
%% occilate too much on the early stopping set
options.stepSize = 2e-07;

%% Use 80% of the data for the training set
trainingIdx = [1:floor(.8*globDat.nSample)];

%% Use the remaining 20% for the stopping set
stoppingIdx = [floor(.8*globDat.nSample)+1:globDat.nSample];

%% Train and viusalize the STRF
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);
preprocRTAC_vis(strfTrained);

%% Well, now the 2nd RTC filters near the peak of the filter time course
%% look pretty decent. But is there any way we can do better? One thing 
%% to consider is that we are using two distinct continuous segments of  
%% a natural movie for the training and stopping sets. The stimulus 
%% correlations within a continuous segment are likely to be quite high.
%% The effect of this is that we are training on a limited stimulus set
%% and testing the stopping criteria against another limited stimulus
%% set. As a potentially better alternative, we could use random samples 
%% from the whole stimulus set as our training set, and use the remaineder
%% for our stopping set. Let's see what happens...

%% create index of all stimuli
idx = [1:globDat.nSample];

%% Create a training index of random samples of 80% of the stimuli
trnIdx = randsample(size(idx,2), .8*size(idx,2));
trainingIdx = idx(trnIdx);

%% Make the stopping index the remaining samples
sIdx = setdiff([1:size(idx,2)],trnIdx);
stoppingIdx = idx(sIdx);

%% To speed this up I have selected a larger step size that works
options.stepSize = 2e-07;

%% Train and viusalize the STRF
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);
preprocRTAC_vis(strfTrained);

%% Not too shabby. This is clearly a better approach. You might be wondering
%% what would happen were you to use different random samples, repeat the
%% sampling a few times, etc. Well that's a great idea, and STRFlab includes
%% functions that take care of all the details of doing multiple resamplings
%% and fittings with your data, meaning you don't have to worry about the 
%% complicated indexing above. First we will look at resampCrossVal, the
%% built-in cross-validation routine

%% First set your options to the defaults for resampCrossVal
options = resampCrossVal;

%% Since the resampling routine now sits between the optimization function
%% and the data, we have to tell resampCrossVal what optimization function
%% to use, in this case trnGradDesc
options.optimOpt = trnGradDesc;

%% set all the options for trnGradDesc
options.optimOpt.earlyStop = 1;
options.optimOpt.nDispTruncate = 0;
options.optimOpt.display = -1;
options.optimOpt.stepSize = 2e-07;

%% now since resampCrossVal will take care of all the indexing complexities
%% we just need to pass an index of all the samples to be used in the
%% resamplings
trainingIdx = [1:globDat.nSample];

%% Now just train the strf. By default this will perform 5-fold cross validation
%% with different training and stopping sets each time. This may take a while...
strfTrained=strfOpt(strf,trainingIdx,options);

%% strfTrained is now a structure containing 5 strf structures within it. You
%% can visualize each one separately
preprocRTAC_vis(strfTrained(1));
preprocRTAC_vis(strfTrained(2)); %% etc...

%% or you can make a new average strf from the mean of all the strf parameters
strf2 = strfTrained(1);
strf2.b1 = mean(cat(2, strfTrained.b1), 2);
strf2.w1 = mean(cat(4, strfTrained.w1), 4);
preprocRTAC_vis(strf2);

%% We see that the 1st RTC filters near the peak of the time course look pretty
%% decent now. It seems we can fit a second-order model (which is a subset of
%% the full RTC as determined by the locality parameter) with natural movies
%% without overfitting to noise and stimulus correlations. This isn't possible 
%% with standard RTC, so again the advantages of STRFlab are obvious. But as 
%% cool as this is, there are limitations of this model that we will consider 
%% in a moment. But first...

%% You will usually want to see how well your STRF model can generalize to
%% new data sets. There are a couple ways to do this. You can initally set the 
%% training index to a large subset of the data you have globalized using
%% strfData, train the STRF model, then use strfFwd with the trained STRF
%% and an index of the remaining data. strfFwd calculates the strf model's
%% output for the indexed stimuli. You can then compare the predicted response
%% to the actual response by taking their correlation and plotting them
%% against each other:
% trainingIdx = [1:floor(.9*globDat.nSample)];
% testingIdx = [floor(.9*globDat.nSample)+1:globDat.nSample];
% strfTrained=strfOpt(strf,trainingIdx,options);
% [strf,predresp]=strfFwd(strfTrained,testingIdx);
% corr(predresp,resp(testingIdx))
% figure; plot(predresp); hold on; plot(resp(testingIdx), 'r')

%% Or you can train on all the globalized data as we were doing earlier
%% and have a separate set of data you load in later to test the
%% STRF model. This may be preferred if you are running up against memory
%% limitations

%% Let's create a separate test data set to test the model fit earlier
%% using the next 5000 frames from the natural movie. We create the fake
%% data just as we did for the trainined data set. If this we real data
%% obviously we would just load in the next set of responses.
rawTestStim = single(mov(1:10,1:10,15001:20000));
rawTestStim = reshape(rawTestStim, [10*10 5000]);
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
respTest = dotdelay(gabor, rawTestStim);
respTest90 = dotdelay(gabor90, rawTestStim);
respTest = sqrt(respTest.^2 + respTest90.^2);
respTest = [zeros([4 1]); respTest(1:end-4)];
respTest(respTest<150) = 0;
respTest(respTest>=300) = 2;
respTest(respTest>=150) = 1;

%% To model a complex cell we turn on RTC in addition to RTA
params.RTAC = [1 1];
params.locality = 1;

%% Preprocess this next set of stimuli
rawTestStim = reshape(rawTestStim, [10 10 5000]);
[stimTest,params] = preprocRTAC(rawTestStim, params);

%% Put the testing stimuli into the global variable so strfFwd can use it
%% to make predicted responses. We put the fake responses in here too, but
%% these are not used because we are not training the model, just looking
%% at the model's output
strfData(stimTest,respTest);

%% Get predicted response from STRF model using strfFwd. We exclude the 
%% first 6 frames as before because our fake responses start on the 7th frame
testingIdx = [1:globDat.nSample];
[strf2,predresp]=strfFwd(strf2,testingIdx);

%% Make sure there are no NaNs in data so we can take a correlation. The first
%% few responses are NaNs because there is a delay between stimulus and response
nonnanidx = find(~isnan(predresp));
corr(predresp(nonnanidx),respTest(nonnanidx))
figure; plot(predresp(nonnanidx)); hold on; plot(respTest(nonnanidx), 'r')

%% We can see that the predictions of our model (in blue) have peaks that line
%% up with most of the spikes in the test data set. The blue line obviously has
%% a different character than the all-or-nothing spikes because it was
%% generated by a linear model rather than a non-linear threshholding


%% So far, to model complex cells we have only been considering covariance
%% across space, that is the covariance terms only refer to pixel interactions
%% within a stimulus frame. This seems to work pretty well around the peak
%% of a 3D filter, but outside of that is noise and we are left with little
%% information about the timecourse of the 3D filter. Previous studies have
%% been limited to considering either covariance in the XY spatial dimensions
%% or in the XT space-time dimensions because of a parameter explosion with
%% each new dimension. By adding the locality parameter, we should be able to
%% cut down the parameters enough to try RTC across the full 3D XYT dimensions.

%% Lets make a complex cell with a preferred direction and a certain temporal
%% frequency
gparams = [0.5000 0.5000 0 5.8654 1.4205 0.0852 0.2816 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);
resp = [zeros([4 1]); resp(1:end-4)];
resp(resp<75) = 0;
resp(resp>=150) = 2;
resp(resp>=75) = 1;

%% To model a complex cell we turn on RTC in addition to RTA
params = preprocRTAC;
params.RTAC = [1 1];

%% Lets leave the locality paramter at 1 for now. This will consider
%% interactions between pixels either adjacent in space or time.
%% Increasing the locality parameter can lead to high memory usage as
%% the parameter number increases
params.locality = 1;

%% Since we are now considering covariance across time as well, we set
%% covtime to 1
params.covtime = 1;

%% This parameter must be set to determine how many frames back to 
%% compute interaction terms for. When covariance across time is not
%% computed this defaults to 0. 
params.covdelays = 7;

%% Make sure rawStim is the proper size before sending to preprocessing
rawStim = reshape(rawStim, [10 10 15000]);
[stim,params] = preprocRTAC(rawStim, params);
strfData(stim,resp)

%% Create new strf. Because the covariance aross time is now considered
%% in the preprocessing, it would be redundant to include delays in the
%% strf, therefore the delay is set to [0]
strf = linInit(size(stim,2), [0]);
strf.b1 = mean(resp);
strf.params = params;

%% Train and visualize the model
options=trnSCG;
options.display=-5;
options.maxIter = 1500;
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);
preprocRTAC_vis(strfTrained);

%% We can see that the 3D filters are recovered for a complex cell. This
%% constrained RTC method is clearly capable of modeling complex cells
%% defined by the Adelson/Bergen energy model in 3D without a huge
%% explosion in the size of the problem.

%% By making some reasonable assumptions and using iterative gradient
%% descent methods, we have pushed RTA and RTC about as far as possible.
%% Other, more complex models may be desireable, however, especially when
%% attempting to model areas beyond V1. Currently our best model consists
%% of preprocessing the stimuli with a pyramid of 3D Gabor filters, that
%% vary in scale, spatial frequency, temporal frequency, orientation and 
%% position which are then either rectified or combined and squared, to give
%% an assortment non-linear filters resembling of simple  and complex cells.
%% The fitting routine determines the linear combination of these filters 
%% that best fits the data. There are many free parameters in this model
%% that determine all aspects of the set of Gabor filters and the 
%% non-linearity. The documentation in preprocWavelets3d explains what 
%% each of these free parameters do.

%% Set the params to the defaults for preprocWavelets3d. Modify and 
%% experiment to see how the size of stim changes. Each column of stim
%% represents the response of a filter to the stimulus
params = preprocWavelets3d;
[stim,params] = preprocWavelets3d(rawStim, params);

%% We'll use the resp derived from the 3D filter in the previous example
strfData(stim,resp)

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% Let's use 5-fold cross validation and coordinate descent with early
%% stopping to regularlize our strf
options = resampCrossVal;
options.optimOpt = trnGradDesc;
options.optimOpt.coorDesc = 1;
options.optimOpt.earlyStop = 1;
options.optimOpt.nDispTruncate = 0;
options.optimOpt.display = -1;
options.optimOpt.stepSize = 0.005;

%% Train strf and average parameters across cross validation runs
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);
strf2 = strfTrained(1);
strf2.b1 = mean(cat(2, strfTrained.b1), 2);
strf2.w1 = mean(cat(4, strfTrained.w1), 4);
preprocWavelets3dVis(strf2);

rawTestStim = reshape(rawTestStim, [10 10 5000]);
[stimTest,params] = preprocWavelets3d(rawTestStim, params);
strfData(stimTest,respTest);

%% Get predicted response from STRF model using strfFwd. We exclude the 
%% first 6 frames as before because our fake responses start on the 7th frame
testingIdx = [1:globDat.nSample];
[strf2,predresp]=strfFwd(strf2,testingIdx);

%% Make sure there are no NaNs in data so we can take a correlation. The first
%% few responses are NaNs because there is a delay between stimulus and response
nonnanidx = find(~isnan(predresp));
corr(predresp(nonnanidx),respTest(nonnanidx))
figure; plot(predresp(nonnanidx)); hold on; plot(respTest(nonnanidx), 'r')


%% Clearly this method can recover the 3D filter. It also has some 
%% significant advantages over the RTC method used above. The number of 
%% parameters needed for fitting in the RTC method will greatly increase 
%% as resolution of the stimulus increases. With the wavelet preprocessing, 
%% however, the number of parameters used in fitting can be controlled by
%% means of the free parameters. Using the wavelet preprocessing, or some
%% other linearizing transform also allows the investigation of a wider
%% range of non-linear transforms of the stimulus as potential models
%% than is possible with RTC


%%% LARS section
options = resampCrossVal;
options.optimOpt = trnLARS;
options.testFrac = 0;
options.nFold = 3;
options.optimOpt.earlyStop=0;
options.optimOpt.maxIter=200;

trainingIdx = [1:globDat.nSample];
[strfTrained, trainoptions, cvResult] =strfOpt(strf,trainingIdx,options);

[preds, maxidx] = optimizeLARSstep(strfTrained, trainoptions, cvResult);

options = trnLARS;
options.earlyStop = 0;
options.maxIter = median(maxidx);

strfTrained=strfOpt(strf,trainingIdx,options);
preprocWavelets3dVis(strfTrained);

strfData(stimTest,respTest);

%% Get predicted response from STRF model using strfFwd. We exclude the 
%% first 6 frames as before because our fake responses start on the 7th frame
testingIdx = [1:globDat.nSample];
[strfTrained,predresp]=strfFwd(strfTrained,testingIdx);

%% Make sure there are no NaNs in data so we can take a correlation. The first
%% few responses are NaNs because there is a delay between stimulus and response
nonnanidx = find(~isnan(predresp));
corr(predresp(nonnanidx),respTest(nonnanidx))
figure; plot(predresp(nonnanidx)); hold on; plot(respTest(nonnanidx), 'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Auditory example

%% specify where wav files and spike files are
base_dir = '/auto/k2/share/strflabGOLD/fakedata/auditory/';

%% use default parameters for the short-time Fourier transform preprocessing
params = preprocSTFT;

%% concatinate the preprocessed stimuli and responses from each stimulus into matrices and
%% keep track of indicies of each stimulus using the assign vector
allstim = [];
allspike = [];
assign = [];
for ii = 1:20
    ii
    [rawStim, fs] = wavread([base_dir, 'stim', num2str(ii), '.wav']);
    params.fs = fs;
    [stim, params] = preprocSTFT(rawStim, params);
    allstim = [allstim; stim];
    
     [spike_time, num_trials] = preprocSpikeTime([base_dir, 'spike', num2str(ii)],size(stim,1));
     allspike = [allspike spike_time];
     assign = [assign ii*ones([1 size(stim,1)])];
end

%% take mean spike rate across all trials
allspike = mean(allspike,1);

global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.

%% Put stimulus and response and group assignments into globDat
strfData(allstim,allspike,assign);

%% Exclude data from the 12th stimulus from training
datIdx =find(assign ~=12);

%% Set options for gradient descent
options=trnGradDesc;
options.display=-1;
options.coorDesc=0;
options.earlyStop=1;
options.stepSize = .00005;
options.nDispTruncate = 0;

%% Initialize strf with up to 40 delays because we are using very small time bins
strf=linInit(size(allstim,2),[0:40]);

%% set the bias term to the mean firing rate
strf.b1 = mean(globDat.resp(datIdx));

%% Use 80% of the data for training
trainingIdx = [1:floor(.8*size(datIdx,2))];

%% Use the remaining 20% for the stopping set
stoppingIdx = [floor(.8*size(datIdx,2))+1:size(datIdx,2)];

%% train the strf
[strfTrained,options]=strfOpt(strf,datIdx(trainingIdx),options, datIdx(stoppingIdx));

%% try prediction using the held out 12th stimulus
datIdxPred =find(assign ==12);
[strfTrained,predResp]=strfFwd(strfTrained,datIdxPred);
sp = allspike(datIdxPred);
corr(sp',predResp)
figure;
plot(sp'); hold on; plot(predResp, 'r')

%% the weights are basically a spectrogram so lets look at them
figure; imagesc(squeeze(strfTrained.w1))


%% Let's trying fitting this STRF with a new (and somewhat experimental) fitting routine
%% trnPF is a algorith that lets you choose an arbitrary Tikhonov regularization matrix, and 
%% decreases the strength of the prior until the STRF's ability to generalize to a stopping
%% set gets worse
options = trnPF;
options.display=-1;
options.maxIter = 1000;
options.lamdaInit = 20000000;

%% We're going to specify a smooth prior the size of our weights so that weights change
%% smoothly in time and space
options.A = full(getSmoothnessPrior([size(squeeze(strf.w1)) 1], [1 1 0]));

%% the last term in the matrix would be for the bias term, but the bias should not be smooth
%% with respect to the other parameters so we set this to 0
options.A(end+1,end+1) = 0;

%% Train and viusalize the STRF
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);


%% try prediction using the held out 12th stimulus
datIdxPred =find(assign ==12);
[strfTrained,predResp]=strfFwd(strfTrained,datIdxPred);
sp = allspike(datIdxPred);
corr(sp',predResp)
figure;
plot(sp'); hold on; plot(predResp, 'r')

%% the weights are basically a spectrogram so lets look at them
figure; imagesc(squeeze(strfTrained.w1))
