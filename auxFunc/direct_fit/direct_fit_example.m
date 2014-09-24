function direct_fit_example(dataDir)

    stimFiles = get_filenames(dataDir, 'stim[0-9]*.wav', 1);
    respFiles = get_filenames(dataDir, 'spike[0-9]*', 1);
    
    timeVaryingPsth = 1;
    dfParams = struct;
    dfParams.timevary_PSTH = timeVaryingPsth;
    
    %default arguments, use STFT
    srData = preprocess_sound(stimFiles, respFiles);
    
    %run direct fit
    strfs = run_direct_fit(srData, dfParams);

    bestStrf = strfs{length(strfs)-1};
    sindx = round(size(bestStrf, 2) / 2);
    bestStrf = bestStrf(:, sindx:size(bestStrf, 2));
    
    %compute predictions for each strf
    modelResponses = compute_predictions(srData, bestStrf, 1, timeVaryingPsth);
    
    for k = 1:length(srData.datasets)
        
        ds = srData.datasets{k}; 
        resp = ds.resp.psth;
        mresp = modelResponses{k};
        
        figure; hold on;
        subplot(2, 1, 1);
        plot_tfrep(ds.stim.tfrep);
        
        subplot(2, 1, 2); hold on;
        plot(resp, 'k-');
        plot(mresp, 'r-');
        legend('Real', 'Model');
        axis tight;
        
    end
    
    %compute coherence for dataset
    [oddPsths, evenPsths] = compute_psth_halves(srData);
    [freq, singleTrialData, psthData, notnormedData] = validate_coherence(oddPsths, evenPsths, modelResponses, srData.respSampleRate);
    
    figure; hold on;
    
	freqCutoff = 50;
	freqint = max(freq) / length(freq);
	freqEnd = round(freqCutoff / freqint);
    indx = 1:freqEnd;
    
    plot(freq(indx), singleTrialData.coherence(indx), 'r-', 'LineWidth', 2);
	plot(freq(indx), psthData.coherence(indx), 'k-', 'LineWidth', 2);
	legend('Model', 'Upper Bound');
	ftitle = sprintf(['coherence: info=%f, upper=%f, lower=%f'], singleTrialData.info, singleTrialData.infoUpper, singleTrialData.infoLower);
	title(ftitle);
	axis([min(freq) freqEnd 0 1]);
    
    
    
    
    
    