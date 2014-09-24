%% Direct Fit training routine
%
%   316/2011. FET Added a new loop to look for STRFs within the lower bound
%   of the maxinformation that have higher tolerance and higher sd.
%   3/16/2011. FET Added option to return all strfs and info values.
%   
%
%   Input:
%       modelParams: 
%
%       trainingIndex: time indicies of data to use as training data,
%           everything else will be held out.
%
%       options:
%           .sparsenesses: a vector of sparseness values to test against
%           .tolerances: a vector of tolerance values to use
%           .separable: 0=space-time separable, 1=non-separable (defaults to 1)
%           .timeVaryingPSTH: whether the time-varying PSTH should be
%               subtracted from each sample (defaults to 0)
%           .timeVaryingPSTHTau: width of hanning window used to smooth
%               time-varying PSTH in ms (defaults to 41ms)
%           .stimSampleRate: post-preprocessing stimulus sampling rate in
%               Hz (defaults to 1000Hz)
%           .respSampleRate: response sampling rate in Hz (defaults to
%               1000Hz. In strflab the sampling rate for the stimulus and
%               the response must be equal.
%           .outputDir: directory where all the stuff goes (defaults to
%               tempdir())
%           .infoFreqCutoff: coherence is computed to judge STRF goodness,
%               this is the cutoff frequency in Hz the coherence is
%               computed to (defaults to 100Hz)
%           .infoWindowSize: size in seconds of the window used to segment
%               responses when computing coherence (defaults to 0.500s)
%           .biasForHighTol: bias the choice of the best STRF towards
%               those with higher tolerance values.
%
%   Output:
%       modelParams: The model Params structure now includes the best STRF 
%       in w1 and in b0. Note that the parameters for w1 are calculated at 
%       all the integer values of from the minimum to the maximum delays given 
%       in modelParams.delays structure.  
%
%       options:
%
function [modelParams, options] = trnDirectFit(modelParams, trainingIndex, options, varargin)
    
    %% set default parameters and return if no arguments are passed
    if nargin == 0
       options = struct;
       options.funcName = 'trnDirectFit';
       options.tolerances = [0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
       options.sparsenesses = [0 1 2 6];
       options.separable = 0;
       options.timeVaryingPSTH = 0;
       options.timeVaryingPSTHTau = 41;
       options.stimSampleRate = 1000;
       options.respSampleRate = 1000;
       options.infoFreqCutoff = 100;
       options.infoWindowSize = 0.500;
       options.storeAllResults = 0;
       options.biasForHighTol = 1;
       
       tempDir = tempname();       
       options.outputDir = tempDir;
       
       modelParams = options;
       
       return;
    end
    
    if ~strcmp(modelParams.type, 'lin')
        error('trnDirectFit only works for linear models!\n');
    end
    
    if (options.respSampleRate ~= options.stimSampleRate)
        error('trnDirectFit: Stimulus and response sampling rate must be equal!\n'); 
    end
    
    global globDat;
    
    fprintf('Writing temp direct fit output to %s\n', options.outputDir);
    [s,mess,messid] = mkdir(options.outputDir);
    
    %% convert strflab's stim/response data format to direct fit's data format
    DS = strflab2DS(globDat.stim, globDat.resp, globDat.groupIdx, options.outputDir, trainingIndex);
    
    %% set up direct fit parameters
    params = struct;

    params.DS = DS;
    params.NBAND = size(globDat.stim, 2);
    params.Tol_val = options.tolerances;
    params.setSep = options.separable;
    params.TimeLagUnit = 'frame';
    %params.timevary_PSTH = options.timeVaryingPSTH;
    params.timevary_PSTH = 0;
    params.smooth_rt = options.timeVaryingPSTHTau;
    params.ampsamprate = options.stimSampleRate;
    params.respsamprate = options.respSampleRate;
    params.outputPath = options.outputDir;
    params.use_alien_space = 0;
    params.alien_space_file = '';

    params.TimeLag = ceil(max(abs(modelParams.delays)));
    
    %% run direct fit
    strfFiles = direct_fit(params);
        
    %% get computed stim and response means 
    svars = load(fullfile(options.outputDir, 'stim_avg.mat'));
    stimAvg = svars.stim_avg;
    respAvg = svars.constmeanrate;
    tvRespAvg = svars.Avg_psth;
    clear svars;
    
    numSamples = length(DS);
    
    %% compute some indicies to use later
    halfIndx = params.TimeLag + 1;   % This is the point corresponding to zero    
    startIndx = halfIndx + round(min(modelParams.delays));
    %endIndx = startIndx + length(modelParams.delays);
    endIndx = halfIndx + round(max(modelParams.delays));
    strfRng = startIndx:endIndx;    
    
    %% subtract mean off of stimulus (because direct fit does this)
    for k = 1:size(globDat.stim, 1)       
        globDat.stim(k, :) = globDat.stim(k, :) - stimAvg';
    end

    %% reorder stims, responses, and group index to accomodate the training index
    allstim = globDat.stim(trainingIndex, :);
    allresp = globDat.resp(trainingIndex);
    groupindx = globDat.groupIdx(trainingIndex);

    groups = unique(groupindx);
    ngindx = zeros(1, length(groupindx));
    for k = 1:length(groups)
       g = groups(k);
       indx = find(groupindx == g);
       ngindx(indx) = k;
    end

    
    %% compute information values for each set of jacknifed strfs per tolerance value
    fprintf('Finding best STRF by computing info values across sparseness and tolerance values...\n');
    
    bestInfoVal = -1;
    bestLowerInfoVal = -1;
    bestStrf = -1;
    bestTol = -1;
    bestSparseness = -1;
    
    % Save each STRF and each info value
    nTol = length(strfFiles);
    nSparse = length(options.sparsenesses);
    nJN = numSamples;
    infos = cell(nTol,nSparse);
    strfs = cell(nTol,nSparse);
    spvals = options.sparsenesses;
    
    for k = 1:length(strfFiles)    %for each tolerance value
        svars = load(strfFiles{k});
        
        %strfsJN is an MxPxT matrix, where M=# of channels, P=# of STRF
        % delays, T = # of stim/response pairs, each element strfsJN(:, :, k) is a STRF
        % constructed from fitting all but pair k to the data.
        strfsJN = svars.STRFJN_Cell;
        
        %strfsJN_std is also an MxPxT matrix. Element strfsJN_std(:, :, k)
        % is the STRF standard deviation across the set of all jacknifed
        % STRFs excluding the kth STRF
        strfsJN_std = svars.STRFJNstd_Cell;
        
        %strfMean is the mean across all jacknifed STRFs for a given
        % tolerance
        strfMean = svars.STRF_Cell;        
        
        %strfStdMean is the average standard deviation across all jacknifed
        % STRFs for a given tolerance
        % strfStdMean = squeeze(mean(strfsJN_std,3));
        strfStdMean = strfsJN_std;
        
        clear svars;
        
        for q = 1:length(spvals)    %for each sparseness value

            %smooth the strf by masking it with a sigmoid-like mask,
            % scaled by the # of deviations specified by the sparseness
            % parameter
            smoothedMeanStrf = df_fast_filter_filter(strfMean, strfStdMean, spvals(q));            
            smoothedMeanStrfToUse = smoothedMeanStrf(:, strfRng);
            
            % Save smoothed STRF
            strfs{k,q} = smoothedMeanStrfToUse;
            
            %the following loop goes through each jacknifed strf for the given tolerance
            %value, smooths it with the sparseness parameter, and then
            %predicts response to it's corresponding held-out stimulus.
            %the coherence and information values are computed and recorded,
            %then the average info value is used to judge the goodness-of-fit
            %for smoothedMeanStrfToUse            

            infoSum = 0; %this will be used to compute the average 
            numJNStrfs = numSamples;
            
            modelResponse = [];
            realResponse = [];
            for p = 1:numJNStrfs
                %jacknifed STRF p (strfsJN(:, :, p)) was constructed by holding out stim/resp pair p
                smoothedMeanStrfJN = df_fast_filter_filter(strfsJN(:, :, p), strfsJN_std, spvals(q));
                strfToUse = smoothedMeanStrfJN(:, strfRng);
		
                %get the held-out stimulus
                srRange = find(ngindx == p);
                stim = allstim(srRange, :);
                rresp = allresp(srRange);                
                gindx = ones(1, size(stim, 1));
		
                %compute the prediction for the held out stimulus
                mresp = conv_strf(stim, modelParams.delays, strfToUse, gindx);

                %add the mean back to the PSTH if necessary                
                if ~options.timeVaryingPSTH
                    mresp = mresp + respAvg;
                else 
                    mresp = mresp + tvRespAvg(p, 1:length(mresp));
                end
                
                modelResponse = cat(1, modelResponse, rv(mresp));
                realResponse = cat(1, realResponse, rv(rresp));
                
            end

            cStruct = compute_coherence_mean(modelResponse, realResponse, options.respSampleRate, options.infoFreqCutoff, options.infoWindowSize);
            avgInfo = cStruct.info;
            infos{k,q} = cStruct;
            
            fprintf('Tolerance=%f, Sparseness=%d, Avg. Prediction Info=%f\n', options.tolerances(k), spvals(q), avgInfo);
            
            %did this sparseness do better?
            if avgInfo > bestInfoVal
            
                bestTol = options.tolerances(k);
                bestSparseness = spvals(q);
                bestInfoVal = avgInfo;
                bestStrf = smoothedMeanStrfToUse;
                bestLowerInfoVal = cStruct.infoLower;
                
            end
            
            clear concatModelResp;
        end
        clear strfsJN;
        clear strfsJN_std;
        clear strfMean;
    end
    
    if options.biasForHighTol
        %% Put a bias for higher tolerance and greater sp
        % Note that bestLowerInfoVal is not updated in the loop. We are looking
        % for any combiniation of tolerance and sparseness that are within the
        % info peak    
        for k = 1:length(strfFiles)                
            for q = 1:length(spvals)
                cStruct = infos{k,q};
                if (cStruct.info > bestLowerInfoVal && options.tolerances(k) >= bestTol && spvals(q) >= bestSparseness)
                    bestTol = options.tolerances(k);
                    bestSparseness = spvals(q);
                    bestInfoVal = cStruct.info;
                    bestStrf = strfs{k,q};
                end
            end
        end
    else
        fprintf('Not biasing for high tolerance...\n');
    end
    
    %% get best strf
    fprintf('Best STRF found at tol=%f, sparseness=%d, info=%f bits\n', bestTol, bestSparseness, bestInfoVal);
    
    modelParams.w1 = bestStrf;
    modelParams.b1 = respAvg;
    modelParams.dfTol = bestTol;
    modelParams.dfSp = bestSparseness;
    
    if (options.storeAllResults)  % Return all results in the internal structure
        modelParams.internal.strfs = strfs;
        modelParams.internal.info = infos;
    end
    
    fprintf('Deleting temp direct fit directory %s...\n', options.outputDir);
    [stat, mess, id] = rmdir(options.outputDir, 's');

    