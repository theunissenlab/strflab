%% Preprocess audio waveforms into strflab-ready time-frequency representation
%
%  Input:
%
%       audioWaveforms: A cell array of digitized sound waveforms or file
%           paths to .wav files.
%
%       params: Structure of preprocessing parameters:
%
%           .rawSampleRate: sample rates of audio waveforms in Hz, not
%               required if audioWaveforms specifies .wav files.
%
%           .tfType: time-frequency type, either 'stft' for short-time
%               Fourier transform, 'lyons' for Lyon's model, or 'wavelet'
%               for wavelet representation.
%
%           .tfParams: time-frequency parameter structure, specific to the tfType
%               specified. For information on these params, see ./sound/timefreq.m.
%               Defaults to default tf-representation params in timefreq.m.
%
%           .outputDir: directory to write all preprocessed stim to, defaults to
%               system temp directory.
%
%           .outputPattern: a pattern to name preprocessed stimuli output
%               files. Must contain '%d' to match the # of the stimulus, for
%               example 'mystim_%d.mat'. Defaults to 'preprocessed_stim_%d.mat'.
%
%           .overwrite: if a cached file is found, ignore it and overwrite it
%
%   Output:
%   
%       wholeStim: an NxM matrix representing the concatenated stimuli,
%           where N is the # of time points and M is the # of features.
%
%       groupIndex: an 1xN vector representing the group index of each
%           time point, i.e. the stimulus/response pair it belongs to.
%
%       stimInfo: information about the preprocessed stimulus:
%           .stimLengths = vector of lengths of each stimulus, in seconds
%           .sampleRate = sample rate of preprocessed stimulus in Hz
%           .numStimFeatures = number of stimulus features (M)
%           .tfType = time frequency representation, same as params.tfType
%           .tfParams = same as params.tfParams;
%           .f = vector of frequencies, one for each feature
%
%       params: same params as passed into function
%
function [wholeStim, groupIndex, stimInfo, params] = preprocSound(audioWaveforms, params)


    %% set default parameters
    if ~isfield(params, 'tfType')
        params.tfType = 'stft';
    end
    
    if ~isfield(params, 'tfParams')
        params.tfParams = struct;
    end
    
    if ~isfield(params, 'outputDir')       
        params.outputDir = tempdir();
    end
    
    if ~isfield(params, 'outputPattern')       
        params.outputPattern = 'preprocessed_stim_%d.mat';
    end
    
    if ~isfield(params.tfParams, 'log')       
        params.tfParams.log = 1;
    end
    
    if ~isfield(params.tfParams, 'dbnoise')       
        params.tfParams.dbnoise = 80;
    end
    
    if ~isfield(params.tfParams, 'refpow')       
        params.tfParams.refpow = 0;
    end
    
    if ~isfield(params, 'overwrite')
        params.overwrite = 0;
    end
    
    if ~isfield(params, 'cache')
        params.cache = 1;
    end
    
    %% check parameters
    allowedTypes = {'stft', 'wavelet', 'lyons'};    
    if ~ismember(params.tfType, allowedTypes)
        error('Unknown time-frequency representation type: %s\n', params.tfType);
    end
    

    %% read .wav files if they're specified instead of audio waveforms
    for k = 1:length(audioWaveforms)        
        if ischar(audioWaveforms{k})
            [adata, sRate, depth] = wavread(audioWaveforms{k});            
            audioWaveforms{k} = adata;
            params.rawSampleRate = sRate;
        end        
    end

    
    %% turn audio waveforms into spectrograms
    f = -1;
    totalStimLength = 0;
    numStimFeatures = -1;
    stimSampleRate = 1000;
    stimStructs = cell(length(audioWaveforms), 1);
    maxPower = -1;
    for k = 1:length(audioWaveforms)    
        
        stim = struct;        
        stimOutputFname = fullfile(params.outputDir, sprintf(params.outputPattern, k));
        fid = fopen(stimOutputFname);
        fileExists = fid ~= -1;
        if fileExists
            fclose(fid);
        end
        if fileExists && ~params.overwrite && params.cache
            %fprintf('Using cached preprocessed stimulus from %s\n', stimOutputFname);
            fvars = load(stimOutputFname);
            stim = fvars.stim;
            clear fvars;
        else
            stim.tfrep = timefreq(audioWaveforms{k}, params.rawSampleRate, params.tfType, params.tfParams);                        
            stim.stimLength = size(stim.tfrep.spec, 2) / stimSampleRate;
            stim.sampleRate = stimSampleRate;
            if params.cache
               save(stimOutputFname, 'stim');
            end
        end        
        stimStructs{k} = stim;
        
        if numStimFeatures == -1
            numStimFeatures = size(stim.tfrep.spec, 1);
        end
        if f == -1
            f = stim.tfrep.f;
        end
        totalStimLength = totalStimLength + size(stim.tfrep.spec, 2);
    
        maxPower = max([maxPower max(stim.tfrep.spec(:))]);        
        
    end
    
    
    %% normalize spectrograms and take log if requested
    if params.tfParams.log        
        if params.tfParams.refpow == 0            
            refpow = maxPower;
        else
            refpow = params.tfParams.refPow;
        end
                
        for k = 1:length(stimStructs)           
            stim = stimStructs{k};
            stim.tfrep.spec = max(0, 20*log10(stim.tfrep.spec/refpow)+params.tfParams.dbnoise);
            stimStructs{k} = stim;
        end        
    end
    
    
    %% concatenate stims into big matrix and record info into struct
    stimInfo = struct;
    stimInfo.stimLengths = zeros(1, length(audioWaveforms));
    stimInfo.sampleRate = stimSampleRate;
    stimInfo.numStimFeatures = numStimFeatures;
    stimInfo.tfType = params.tfType;
    stimInfo.tfParams = params.tfParams;
    stimInfo.f = f;
    
    groupIndex = zeros(1, totalStimLength);
    wholeStim = zeros(totalStimLength, numStimFeatures);
    
    cindx = 1;
    for k = 1:length(stimStructs)
       
        stim = stimStructs{k};
        slen = size(stim.tfrep.spec, 2);
        tend = cindx + slen - 1;
        trng = cindx:tend;
        
        wholeStim(trng, :) = stim.tfrep.spec';
        groupIndex(trng) = k;
        stimInfo.stimLengths(k) = slen / stimSampleRate;
        
        cindx = tend + 1;
        
    end
    
    clear stimStructs;
    