%% make some test data
stim=randn(20,20,3000);

h=make3dgabor([20 20 1], [0.5 0.5 90 2 0 0.1 0.1 0]);

h = randn(20,20);

for ii=1:3000
  resp(ii,1)=sum(sum(h.*stim(:,:,ii)))^2+.5*randn;
end
resp = resp/max(resp(:))*2;
resp = [0 0 resp(1:end-2)']';  % some time shifts

resp(100:120) = NaN;  % add NaNs
resp(300:310) = NaN;


%% physical value settings (optional)
physics.binwidth_t = 1000/83;
physics.binwidth_s = 1/18 * 5;
physics.unit_t = 'ms';
physics.unit_s = 'deg';
physics.centerx = -10;
physics.centery = 5;


%% preprocessing settings
params = preprocWavelets;
params.phasemode = 2;  % half-rectified sin and cos amplitudes
params.sfmax = 5; % maximum spatial freq.: 5 cycles/ stimulus size
params.std_step = 3.0;
  % separation between the adjacent wavelets in STD of Gaussian window

%% !! Uncomment the following 2 lines to test static Gabor wavelets !!
%  params.veldivisions = 1; %% static Gabor mode
%  params.tsize = 1; %% static Gabor mode


%% preprocess
[stim params] = preprocWavelets(stim, params);


%% model settings
delays = [0:6];
strf = linInit(params.nChan, delays);

strf.params = params;
strf.physics = physics;

strf.w1=strf.w1*0;
strf.b1=nanmean(resp);


%% optimizer settings
options=trnGradDesc;

options.coorDesc=1;
options.display=1;
options.maxIter = 50;
options.earlyStop = 1;

global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(stim,resp);

trainingIdx=1:floor(9*globDat.nSample/10);  % generate index of training samples.
earlyStopIdx=(trainingIdx+1):globDat.nSample;  % generate index of early stopping samples.


%% optimize
[strf options] = strfOpt(strf, trainingIdx, options, earlyStopIdx);


%% visualize
visopt.show_bests = 1;
preprocWaveletsVis(strf, visopt);
