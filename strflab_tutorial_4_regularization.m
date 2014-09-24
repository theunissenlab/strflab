%% Now that we have the STRFLab basics down, let's look more in depth
%% at the topic of regularization and how it can help STRF estimates
%% when the response of the system is noisy

%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_4_regularization');
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

%% Let's create a complex gabor filter and the response as before
gparams = [.5 .5 0 5.5 0 .0909 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 5], gparams);
gabor = reshape(gabor, [10*10 5]);
gabor90 = reshape(gabor90, [10*10 5]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);
resp = sqrt(resp.^2 + resp90.^2);

%% Now let's add noise to the computed response. We will use gaussian 
%% noise with a mean of 0 and a standard deviation 1.5 times the 
%% standard deviation of the original response. We will see how well 
%% we can recover the filters given this large amount of noise
resp = resp + 1.5*std(resp)*randn(size(resp));

%% delay the response as before
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

%% Let's increase the number of pixel covariances included in the model
%% by setting the locality parameter higher. This will, however, make the
%% model more succeptible to fitting parameter weights to noise and thus
%% yielding poor filter estimates. 
params.locality = 2;

%% Make sure rawStim is the proper size before sending to preprocessing
rawStim = reshape(rawStim, [10 10 15000]);
[stim,params] = preprocRTAC(rawStim, params);
strfData(stim,resp)

%% Create new strf
strf = linInit(size(stim,2), [0:8]);
strf.b1 = mean(resp);
strf.params = params;

%% We will use SCG as before to descend the error surface quickly and 
%% set the maximum number of iterations high to make sure it converges
options=trnSCG;
options.display=-5;
options.maxIter = 1500;
trainingIdx = [1:globDat.nSample];
strfTrained=strfOpt(strf,trainingIdx,options);

%% Now lets see what we have as far as filter estimates
preprocRTAC_vis(strfTrained);

%% Generally it looks pretty ugly. The first or second principle 
%% components near the filters peak look vaguly like the filter
%% but are not terribly convincing. Because we are completely descending 
%% the error surface with SCG we are likely fitting to noise in the
%% data. Perhaps we could prevent such overfitting by stopping the
%% decent when it becomes clear we are fitting to noise in the data.
%% This scheme is called early-stopping and can prevent over-
%% fitting to noise in the data and stimulus correlations. Gradient 
%% descent fits the parameters that account for the most variance first,
%% and since most of the variance in natural images is low frequency, 
%% early stopping is effectively low pass filtering the result, making
%% it look better. Can we do early stopping in a principled way rather 
%% than just picking an arbitrary number of iterations for SCG? Yes, but
%% not with SCG. trnSCG tries to descend the gradient as 
%% quickly as possible, taking steps of various sizes and is not well suited 
%% to doing early stopping. So we turn now to another fitting routine 
%% trnGradDesc. This function can take an index of data samples to use as 
%% a stopping set. The gradient will be computed on the training set as 
%% usual and a small step of a specified size will be taken down the 
%% gradient. The error is also assessed at each step on the stopping set 
%% and if the error goes up too much on this separate set, the training 
%% will stop. As shown in the appendix of the paper on STRFlab,
%% early stopping gradient descent corresponds to a Gaussian prior (with 
%% the proper covariance matrix) over the model parameters.

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
options.stepSize = 5e-04;

%% We'll try using 80% of the data for the training set
trainingIdx = [1:floor(.8*globDat.nSample)];

%% And the remaining 20% for the stopping set
stoppingIdx = [floor(.8*globDat.nSample)+1:globDat.nSample];

%% Train and viusalize the STRF
strfTrained=strfOpt(strf,trainingIdx,options,stoppingIdx);

%% Now let's see what early stopping has done for us...
preprocRTAC_vis(strfTrained);

%% The result should look somewhat better with the first two principle
%% components around the filters peak resembling the even and odd phases
%% of a complex gabor filter. It's still noisy to be sure, but a signifcant
%% improvement. Taking advantage of such regularization techniques is a huge
%% advantage the Bayesian regression framework has over the standard way of 
%% doing RTA and RTC. As noted before, the essential change we have introduced
%% is an assumed prior distribution of the parameter weights. Such prior 
%% assumptions can significantly benefit fitting models with noisy data. I 
%% leave it to the user to try adjusting the locality parameter up and down
%% to see the regularization effects of this parameter with and without
%% early stopping.