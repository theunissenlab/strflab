function [surpriseStimLouder, surpriseStimQuieter, groupIndex, stimInfo, params] = preprocSoundSurprise(wholeStim, groupIndex, params)

    %% set default parameters
    if ~isfield(params, 'domainFrequencyWidth')
        params.domainFrequencyWidth = 3;
    end
    
    if ~isfield(params, 'domainTimeWidth')
        params.domainTimeWidth = 3;
    end
    
    if ~isfield(params, 'domainGap')
        params.domainGap = 4;
    end
    
    if ~isfield(params, 'outputPath')
        params.outputPath = tempdir();
    end
    
    if ~isfield(params, 'outputDesc')
        params.outputDesc = 'default';
    end
    
    if ~isfield(params, 'cache')
        params.cache = 0;
    end
    
    
    outputFileName = fullfile(params.outputPath, sprintf('surprise.%s.mat', params.outputDesc));
    if params.cache && exist(outputFileName, 'file')
        vars = load(outputFileName);
        surpriseStimLouder = vars.surpriseStimLouder;
        surpriseStimQuieter = vars.surpriseStimQuieter;
        stimInfo = vars.stimInfo;
        params = vars.params;
        
        return;        
    end
    
    
    %% initialize stuff    
    freqWidth = params.domainFrequencyWidth;
    timeWidth = params.domainTimeWidth;
    timeGap = params.domainGap;

    familiarity = ones(size(groupIndex));
    
    nChannels = size(wholeStim, 2);
    
    cStart = freqWidth + 1;
    cEnd = nChannels - 2*freqWidth;
    
    groups = unique(groupIndex);
    lengths = zeros(1, length(groups));
    for k = 1:length(groups)       
        lengths(k) = sum(groupIndex == k);
    end
    
    use_stim = ones(1, length(groupIndex));
    stim_names = cell(1, length(groups));
    for k = 1:length(groups)
        stim_names{k} = sprintf('stim_%d', k);
    end
    
    surpriseStimLouder = zeros(cEnd-cStart+1, size(wholeStim, 1)); 
    surpriseStimQuieter = zeros(cEnd-cStart+1, size(wholeStim, 1)); 
   
    freqIndex = 1;
    
    %% compute surprise, frequency by frequency
    for freqBand = cStart:cEnd
    
        tic;
        %% Transform into the PC space + target of the stimuli
        try            
            [PCA_Tmat, E, D] = surprise_preprocess(wholeStim', freqBand, timeWidth, timeGap, freqWidth);            
        catch            
            disp('Problem doing PCA on spectrogram; probably out of range.');
            disp(lasterr);
            return
        end
        
        PC_Nbins = [15 4 6 6 6 6 15]; % How many bins to use for S and each kept PC.  The first entry is for S, and the next is for the *weakest* kept PC.
        n_dims_minus_1 = length(PC_Nbins)-2;

        temp_PCA_mat = PCA_Tmat([1 (end-n_dims_minus_1):end],:,1);  % Just the parts of PCA_Tmat that will be used.
        clear PCA_Tmat  %Free up some memory
        PC_Maxes = max(temp_PCA_mat,[],2);
        PC_Mins = min(temp_PCA_mat,[],2);
        
        %% construct the joing density P(S,D)
        [joint_p,which_cubes,what_position] = surprise_get_joint_p_mat(temp_PCA_mat, PC_Maxes, PC_Mins, PC_Nbins, familiarity);
        
        %% construct the conditional density P(S | D)
        cond_p = surprise_get_7D_cond_p(joint_p);
        clear joint_p;
        
        %% convert the conditional density to output 'surprise-o-gram'
        out = faster_surprise_cond_p_7D_to_output(cond_p, which_cubes, what_position, use_stim, lengths, stim_names);        

        sindx = 1;
        for k = 1:length(groups)
            len = lengths(k);
            eindx = sindx + len - 1;
            srng = sindx:eindx;
            surpriseStimLouder(freqIndex, srng) = out(k).louder;
            surpriseStimQuieter(freqIndex, srng) = out(k).quieter;                    
            sindx = eindx + 1;
        end

        freqIndex = freqIndex + 1;        
        etime = toc;
        fprintf('Processed frequency band %d in %0.0fs\n', freqBand, etime);
        
    end
    
    stimInfo = struct;
    
    if params.cache
        save(outputFileName, 'surpriseStimLouder', 'surpriseStimQuieter', 'groupIndex', 'stimInfo', 'params');
    end
    