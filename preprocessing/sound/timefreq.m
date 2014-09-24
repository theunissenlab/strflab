%% General purpose time-frequency representation function
%
%   Input:
%       audioWaveform: digitized sound pressure waveform
%
%       sampleRate: rate in Hz of digitized waveform
%
%       typeName: 'stft' for short-time fourier transforms
%                 'wavelet' for wavelet transforms
%                 'lyons' for lyons-model
%
%       params: depends on typeName, default values used if not given
%
%           'stft' params
%               .fband: spacing between frequency bands of spectrogram (125)
%               .nstd: # of std deviations that define width of Gaussian
%                   window (6)
%               .low_freq: lowest frequency in Hz (250)
%               .high_freq: highest frequency in Hz (8000)
%               .log: take base 10 log of spectrogram (1)
%               .refpow: reference power when taking log spectrogram
%               .dbnoise: noise floor in dB, everything below this is
%                   ignored ( when .log=1)
%
%           'wavelet' params
%               (currently not implemented)
%
%           'lyons' params (Requires AuditoryToolbox)               
%               .low_freq: lowest frequency in Hz (250)
%               .high_freq: highest frequency in Hz (8000)
%               .earQ: quality factor of each filter (8)
%               .agc: use adaptive gain control (1)
%               .differ: use differential gain control (1)
%               .tau: time constant of gain control (3)
%               .step: 1/step is approximately the number of filters per bandwidth 
%
%   Output:
%
%       tfrep: time-frequency representation
%           .type: the name of the type of time frequency representation ('stft', 'lyons', 'wavelet')
%           .t: vector of time points for tf-representation
%           .f: vector of frequencies for tf-representation
%           .spec: matrix of values for tf-representation at a given
%               time-frequency cell
%           .params: the parameters that created the time-frequency representation
%
function tfrep = timefreq(audioWaveform, sampleRate, typeName, params)

    if nargin < 3
        tfrep = make_tfrep(typeName);
    else
        tfrep = make_tfrep(typeName, params);
    end

    tfrep.params.rawSampleRate = sampleRate;
    
    %% create spectrogram
    switch typeName
       
        case 'stft'
            
            %compute raw complex spectrogram
            
            specSampleRate = 1000;
            desiredSampleRate = ceil(sampleRate/specSampleRate)*specSampleRate;
            if desiredSampleRate > sampleRate
                audioWaveform = resample(audioWaveform, desiredSampleRate, sampleRate);
                sampleRate = desiredSampleRate;
            end
            
            twindow = tfrep.params.nstd/(tfrep.params.fband*2.0*pi);   % Window length
            winLength = fix(twindow*sampleRate);  % Window length in number of points
            winLength = fix(winLength/2)*2; % Enforce even window length
            increment = floor(sampleRate/specSampleRate); % Sampling rate of spectrogram in number of points - set at specSampleRate
            
            [s, t0, f0, pg] = GaussianSpectrum(audioWaveform, increment, winLength, sampleRate); 
                       
            %normalize the spectrogram within the specified frequency range
            maxIndx = find(f0 >= tfrep.params.high_freq);
            maxIndx = maxIndx(1);
            minIndx = find(f0 < tfrep.params.low_freq);
            minIndx = minIndx(end) + 1;
            
            normedS = abs(s(minIndx:maxIndx, :));
            
            %set tfrep values
            fstep = f0(2);
            tfrep.t = t0;
            tfrep.f = f0(minIndx):fstep:f0(maxIndx);
            tfrep.spec = normedS;
            
        case 'wavelet'
            fprintf('Wavelets not currently implemented!\n');
            return;
            
        case 'lyons'
            df = sampleRate / 1000; % Decimation factor
            ldata = LyonPassiveEar_new_mod(audioWaveform, sampleRate, df, tfrep.params.low_freq,...
                                           tfrep.params.high_freq, tfrep.params.earQ, tfrep.params.step,...
                                           tfrep.params.differ, tfrep.params.agc, tfrep.params.tau);
            lspec = ldata.spec;
                                       
            tlen = size(lspec, 2);            
            t0 = (0:(tlen-1))*1e-3;
            %finc = round((tfrep.params.high_freq - tfrep.params.low_freq) / flen);
            %f0 = tfrep.params.low_freq:finc:tfrep.params.high_freq;
            f0 = ldata.freqs;
            
            tfrep.spec = flipud(lspec);
            tfrep.t = t0;
            tfrep.f = f0;
    end
    
    