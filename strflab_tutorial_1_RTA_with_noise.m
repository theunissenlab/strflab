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


%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_1_RTA_with_noise');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))

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
