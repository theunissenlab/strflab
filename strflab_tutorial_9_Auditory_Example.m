%% Auditory example: Here we'll do preprocessing of auditory data, and then
%  compare 3 linear fitting methods: Direct Fit, Gradient Descent, and
%  Coordinate Descent

%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_9_Auditory_Example');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))


%% get list of stim files
dataDir = fullfile(strflabDir, 'fakedata', 'auditory');
audioFiles = get_filenames(dataDir, 'stim[0-9]*.wav', 1);


%% preprocess the sound to produce spectrograms
preprocStimParams = struct;      %create preprocessing param structure
preprocStimParams.tfType = 'stft'; %use short-time FT
tfParams = struct;               %create time-frequency params
tfParams.high_freq = 8000;       %specify max freq to analyze
tfParams.low_freq = 250;         %specify min freq to analyze
tfParams.log = 1;                %take log of spectrogram
tfParams.dbnoise = 80;           %cutoff in dB for log spectrogram, ignore anything below this
tfParams.refpow = 0;             %reference power for log spectrogram, set to zero for max of spectrograms across stimuli
preprocStimParams.tfParams = tfParams;

% make a temporary directory to store preprocessed sound files (should be
%  specific to parameters for preprocSound)
tempPreprocDir = tempname();    
[s,mess,messid] = mkdir(tempPreprocDir);
preprocStimParams.outputDir = tempPreprocDir;

[wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(audioFiles, preprocStimParams);


%% preprocess spike times to produce PSTHs
respFiles = get_filenames(dataDir, 'spike[0-9]*', 1);
allSpikeTrials = cell(length(respFiles), 1);
for k = 1:length(respFiles)
    allSpikeTrials{k} = read_spikes_from_file(respFiles{k});
end

preprocRespParams = struct;         %create preprocessing struct
preprocRespParams.units = 'ms';     %our files have units of ms
preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
preprocRespParams.stimLengths = stimInfo.stimLengths; %needed to compute correct PSTH lengths

[wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


%% display stim and responses
displayStimResp = 0;

if displayStimResp
    nGroups = length(unique(groupIndex));
    for k = 1:nGroups

        tRng = find(groupIndex == k);

        stim = wholeStim(tRng, :);
        resp = wholeResponse(tRng);

        tInc = 1 / stimInfo.sampleRate;
        t = 0:tInc:(size(stim, 1)-1)*tInc;

        figure; hold on;
        
        %plot spectrogram
        subplot(2, 1, 1);
        imagesc(t, stimInfo.f, stim'); axis tight;
        axis xy;
        v_axis = axis;
        v_axis(1) = min(t); v_axis(2) = max(t);
        v_axis(3) = min(stimInfo.f); v_axis(4) = max(stimInfo.f);
        axis(v_axis);
        xlabel('Time (s)'); ylabel('Frequency (Hz)');

        %plot PSTH
        subplot(2, 1, 2);
        plot(t, resp, 'k-'); axis([0 max(t) 0 1]);
        xlabel('Time (s)'); ylabel('P[spike]');

        title(sprintf('Pair %d', k));
    end
end


%% Initialize strflab global variables with our stim and responses
global globDat
strfData(wholeStim, wholeResponse, groupIndex);


%% Initialize a linear model that extends 75ms back in time
strfLength = 40;
strfDelays = 0:(strfLength-1);
modelParams = linInit(stimInfo.numStimFeatures, strfDelays);


%% pick training datasets for DirectFit
trainingGroups = 1:18;
trainingIndex = findIdx(trainingGroups, groupIndex);


%% Initialize and run the DirectFit training routine
optOptions = trnDirectFit();
fprintf('\nRunning Direct Fit training...\n');
[modelParamsDF, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);


%% pick training and early stopping datasets for Gradient and Coordinate Descent
trainingGroups = 1:17;
earlyStoppingGroups = [18];

trainingIndex = findIdx(trainingGroups, groupIndex);
earlyStoppingIndex = findIdx(earlyStoppingGroups, groupIndex);


%% Initialize and run Gradient Descent w/ early stopping
optOptions = trnGradDesc();
optOptions.display = 1;
optOptions.maxIter = 1000;
optOptions.stepSize = 2e-6;
optOptions.earlyStop = 1;
optOptions.gradNorm = 1;

%initialize bias term to mean of response
modelParams.b1 = mean(wholeResponse(trainingIndex));

fprintf('\nRunning Gradient Descent training...\n');
[modelParamsGradDesc, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, earlyStoppingIndex);


%% Initialize and run Coordinate Descent w/ early stopping
optOptions = trnGradDesc();
optOptions.display = 1;
optOptions.maxIter = 300;
optOptions.stepSize = 1e-4;
optOptions.earlyStop = 1;
optOptions.coorDesc = 1;

fprintf('\nRunning Coordinate Descent training...\n');
[modelParamsCoorDesc, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, earlyStoppingIndex);


%% split original spike trials into two PSTHs for purposes of validation
preprocRespParams.split = 1;
[wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


%% compute responses to validation data for each STRF
validationGroups = [19 20];
respReal = [];
respRealHalf1 = [];
respRealHalf2 = [];
respDF = [];
respGradDesc = [];
respCoorDesc = [];
for k = 1:length(validationGroups)
   
    g = validationGroups(k);
    gIndx = find(globDat.groupIdx == g);
    stim = globDat.stim(gIndx);
    resp = globDat.resp(gIndx);
    respH1 = wholeSplitResponse(1, gIndx);
    respH2 = wholeSplitResponse(2, gIndx);
        
    [modelParamsDF, resp1] = strfFwd(modelParamsDF, gIndx);
    [modelParamsGradDesc, resp2] = strfFwd(modelParamsGradDesc, gIndx);
    [modelParamsCoorDesc, resp3] = strfFwd(modelParamsCoorDesc, gIndx);
    
    respReal = [respReal resp];
    respRealHalf1 = [respRealHalf1 respH1];
    respRealHalf2 = [respRealHalf2 respH2];
    respDF = [respDF resp1'];
    respGradDesc = [respGradDesc resp2'];
    respCoorDesc = [respCoorDesc resp3'];
end


%% rescale model responses
respDF = (respDF / max(respDF))*max(respReal);
respGradDesc = (respGradDesc / max(respGradDesc))*max(respReal);
respCoorDesc = (respCoorDesc / max(respCoorDesc))*max(respReal);


%% Compute performance on validation datasets for each STRF
avgNumTrials = mean(respInfo.numTrials); %taking mean isn't necessary here
infoFreqCutoff = 90; %Hz
infoWindowSize = 0.500; %500ms
[cBound, cDF] = compute_coherence_full(respDF, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
[cBound, cGradDesc] = compute_coherence_full(respGradDesc, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
[cBound, cCoorDesc] = compute_coherence_full(respCoorDesc, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);


%% plot STRFS
figure; hold on;

subplot(3, 1, 1); hold on;
imagesc(strfDelays, stimInfo.f, modelParamsDF.w1); axis tight;
absmax = max(max(abs(modelParamsDF.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Direct Fit | bias=%f', modelParamsDF.b1));

subplot(3, 1, 2); hold on;
imagesc(strfDelays, stimInfo.f, squeeze(modelParamsGradDesc.w1)); axis tight;
absmax = max(max(abs(modelParamsGradDesc.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Gradient Descent | bias=%f', modelParamsGradDesc.b1));

subplot(3, 1, 3); hold on;
imagesc(strfDelays, stimInfo.f, squeeze(modelParamsCoorDesc.w1)); axis tight;
absmax = max(max(abs(modelParamsCoorDesc.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Coordinate Descent | bias=%f', modelParamsCoorDesc.b1));


%% display predictions
displayPredictions = 1;
if displayPredictions
    
    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(length(respReal)-1)*tInc;
   
    figure; hold on;
    
    subplot(3, 1, 1); hold on;
    plot(t, respReal, 'k-');
    plot(t, respDF, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Direct Fit');
    legend('Real', 'Model');
    
    subplot(3, 1, 2); hold on;
    plot(t, respReal, 'k-');
    plot(t, respGradDesc, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Gradient Descent');
    legend('Real', 'Model');
    
    subplot(3, 1, 3); hold on;
    plot(t, respReal, 'k-');
    plot(t, respCoorDesc, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Coordinate Descent');
    legend('Real', 'Model');
    
end


%% plot coherences and information
figure; hold on;
plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
plot(cDF.f, cDF.c, 'b-');
plot(cGradDesc.f, cGradDesc.c, 'g-');
plot(cCoorDesc.f, cCoorDesc.c, 'r-');
axis tight;
title('Coherences for Validation Set');
legend('Upper Bound', 'Direct Fit', 'Grad Desc', 'Coord Desc');


%% compute performance ratios
perfDF = cDF.info / cBound.info;
perfGradDesc = cGradDesc.info / cBound.info;
perfCoorDesc = cCoorDesc.info / cBound.info;

fprintf('Performance Ratios:\n');
fprintf('\tDirect Fit: %f\n', perfDF);
fprintf('\tGradient Descent: %f\n', perfGradDesc);
fprintf('\tCoordinate Descent: %f\n', perfCoorDesc);
