function strfFiles = direct_fit(params)

global DF_PARAMS
DF_PARAMS = params;

DS = params.DS;
NBAND = params.NBAND;
Tol_val = params.Tol_val;
setSep = params.setSep;
TimeLag = params.TimeLag;
TimeLagUnit = params.TimeLagUnit;
timevary_PSTH = params.timevary_PSTH;
smooth_rt = params.smooth_rt;
ampsamprate = params.ampsamprate;
respsamprate = params.respsamprate;

%the intermediate result path
outputPath = params.outputPath;

% OnlyOne flag used for avoiding overfitting
OnlyOne = 1;

split_data_set = 0;

if (split_data_set)    
     % Divide the orignal data files into two parts

     stim_env = Check_And_Load(DS{1}.stimfiles);
     psth_rec = Check_And_Load(DS{1}.respfiles);
     oldLength = min(length(stim_env), length(psth_rec));

     partI = stim_env(:, 1:floor(0.9*oldLength));
     partII = stim_env(:, floor(0.9*oldLength)+1:oldLength);
     [p, sn, se] = fileparts(DS{1}.stimfiles);
     save(fullfile(outputPath,['partI_' sn,'.mat']), 'partI');
     save(fullfile(outputPath,['partII_' sn, '.mat']), 'partII');

     DS{1}.stimfiles = fullfile(p, ['partI_', sn, '.mat']);
     predDS{1}.stimfiles = fullfile(p, ['partII_', sn, '.mat']);

     DS{1}.nlen = floor(0.9 * oldLength);
     predDS{1}.nlen = oldLength - DS{1}.nlen;
     predDS{1}.ntrials = DS{1}.ntrials;

     % Divide orignal response files

     [p, rn, re] = fileparts(DS{1}.respfiles);
     partI = psth_rec(:, 1:floor(0.9*oldLength));
     partII = psth_rec(:, floor(0.9*oldLength)+1:oldLength);

     save(fullfile(outputPath,['partI_', rn, '.mat']), 'partI');
     save(fullfile(outputPath,['partII_', rn, '.mat']), 'partII');

     DS{1}.respfiles = fullfile(p, ['partI_', rn, '.mat']);
     predDS{1}.respfiles = fullfile(p, ['partII_', rn, '.mat']);
end


% Begin to calculate

if isempty(TimeLagUnit)
    fprintf( 'Please specify Time-lag unit (msec or frames).');
end

if isempty(TimeLag) | isempty(Tol_val) ...
        | isempty(setSep) | isempty(timevary_PSTH)...
        | isempty(smooth_rt)
    fprintf( 'Please set all calculation parameter values.\n');
end

% =========================================================
% calculate avg. of stimulus and response that used later
% =========================================================

[stim_avg, avg_psth, psth, errFlg] = df_cal_AVG(DS);

% Check if cal_Avg ends normally

if errFlg == 1
    fprintf( 'df_cal_AVG ended with error!\n');
end

% =========================================================
% Now calcualting stimulus AutoCorr.
% =========================================================


if isempty(ampsamprate)
    ampsamprate = 1000;
end

if strcmp(TimeLagUnit, 'msec')
    twindow = [-round(TimeLag*ampsamprate/1000) round(TimeLag*ampsamprate/1000)];
elseif strcmp(TimeLagUnit, 'frame')
    twindow = [-TimeLag TimeLag];
end

if setSep == 0  % Nonseparable space-time algorithm
    disp('Now calculating stim autocorrelation');
    do_long_way = 1;
    [cached_dir,maxsize] = df_dir_of_caches;
    autocorr_start_time = cputime;
    hashes_of_stims = df_create_stim_cache_file(outputPath,DS);

    if ~strcmp(cached_dir,'null')
        do_long_way = 0;
        [loaded,order] = sort(hashes_of_stims);  %  Sort to make the checksum invarient to dataset shuffling
        n_trial_array = get_ntrials(DS);
        checksum_for_autocorr_calc = df_checksum(df_load_function_text('df_cal_AutoCorr'),loaded,n_trial_array(order),stim_avg,twindow,NBAND);  % Sort the ntrial array the same way as the stimuli
        [CS,errFlg] = do_cached_calc_checksum_known('df_cal_AutoCorr',checksum_for_autocorr_calc,1, DS, stim_avg, twindow, NBAND);
    else
        [CS, errFlg] = df_cal_AutoCorr(1, DS, stim_avg, twindow, NBAND);
    end

    %[DS_data,the_checksum] = load_DS_data(DS,stim_avg,twindow, NBAND);
    %autocorr_start_time = cputime;
    %[CS, errFlg] = do_cached_calc_checksum_known('df_cal_AutoCorr_for_cache',the_checksum, DS_data, stim_avg,twindow, NBAND);
    autocorr_end_time = cputime;
    disp(['The autocorrelation took ' num2str(autocorr_end_time - autocorr_start_time) ' seconds.']);
    currentPath = pwd;    
    if ~isempty(outputPath)
        cd (outputPath);
    else
        disp('Saving output to Output Dir.');
        stat = mkdir('Output');
        cd('Output');
        outputPath = pwd;
    end

    save('Stim_autocorr.mat', 'CS');
    cd(currentPath);

    % Check if df_cal_AutoCorr ends normally
    if errFlg == 1
        fprintf( 'df_cal_AutoCorr ended in failure!\n');
    end

    % Done calcualtion of stimulus AutoCorr
    % =========================================================
    % Now calcualting stimulus spike Cross Corr.
    % =========================================================

    %  Let's assume that if the user has caching on and is using the GUI
    %  that they might want to evaluate df_cal_CrossCorr more than once; so:
    cache_crosscorr = ~strcmp(cached_dir,'null');  % This is the flag to say if we'll cache results specific to the current spike train.

    hashes_of_stims = df_create_stim_cache_file(outputPath,DS);
    hashes_of_spikes = df_create_spike_cache_file(outputPath,DS);
    if isempty(smooth_rt)
        smooth_rt = 41;
    end
    if ~exist('psth_option')
        if timevary_PSTH == 0
            psth_option = 0;
        else
            psth_option = 1;
        end
    end
    checksum_CrossCorr = df_checksum(df_load_function_text('df_cal_CrossCorr'),hashes_of_spikes,hashes_of_stims,twindow,smooth_rt,psth_option);

    if ~strcmp(cached_dir,'null')
        %[CSR, CSR_JN, errFlg]= do_cached_calc_checksum_known('df_cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        [CSR, CSR_JN, errFlg]= df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        save(fullfile(outputPath,'StimResp_crosscorr.mat'), 'CSR');
        save(fullfile(outputPath,'SR_crosscorrJN.mat'), 'CSR_JN');
    else
        [CSR, CSR_JN, errFlg]= df_cal_CrossCorr(DS,stim_avg,avg_psth,...
            psth,twindow,NBAND);
    end

    % Check if df_cal_CrossCorr ends normally
    if errFlg == 1
        set(handles.figure1, 'Pointer', 'Arrow');
        return
    end
    % =========================================================
    %Done calcualtion of stimulus-spike CrossCorr in GUI window
    % =========================================================
    disp('Calculating strfs for each tol value.');
    %checksum_CrossCorr
    df_calStrf_script;
    calculation_endtime = cputime;
    disp(['The STRF calculation took ' num2str(calculation_endtime - autocorr_start_time) ' seconds.']);

else % Separable space-time algorithm

    % Provide Space-time separability algorithm to estimate STRF

    %         [CSspace, CStime, errFlg] = df_cal_AutoCorrSep(DS, stim_avg,...
    %         twindow, NBAND, 1);


    disp('Now calculating stim autocorrelation');
    do_long_way = 1;
    [cached_dir,maxsize] = df_dir_of_caches;
    autocorr_start_time = cputime;
    hashes_of_stims = df_create_stim_cache_file(outputPath,DS);


    if ~strcmp(cached_dir,'null')
        do_long_way = 0;
        [loaded,order] = sort(hashes_of_stims);  %  Sort to make the checksum invarient to dataset shuffling
        n_trial_array = get_ntrials(DS);
        checksum_for_autocorr_calc = df_checksum(df_load_function_text('df_cal_AutoCorrSep'),'df_cal_AutoCorrSep',loaded,n_trial_array(order),stim_avg,twindow,NBAND);  % Sort the ntrial array the same way as the stimuli
        [CSspace, CStime, errFlg] = do_cached_calc_checksum_known('df_cal_AutoCorrSep',checksum_for_autocorr_calc, DS, stim_avg, twindow, NBAND,1);
    else
        [CSspace, CStime, errFlg] = df_cal_AutoCorrSep(DS, stim_avg, twindow, NBAND,1);
    end

    %[DS_data,the_checksum] = load_DS_data(DS,stim_avg,twindow, NBAND);
    %autocorr_start_time = cputime;
    %[CS, errFlg] = do_cached_calc_checksum_known('df_cal_AutoCorr_for_cache',the_checksum, DS_data, stim_avg,twindow, NBAND);
    autocorr_end_time = cputime;
    disp(['The autocorrelation took ' num2str(autocorr_end_time - autocorr_start_time) ' seconds.']);

    % Check if df_cal_AutoCorrSep ends normally
    if errFlg == 1
        fprintf( 'df_cal_AutoCorrSep ended with error!\n');
    end

    % Calculate cross-correlation between stimuli and spike


    %  Let's assume that if the user has caching on and is using the GUI
    %  that they might want to evaluate df_cal_CrossCorr more than once; so:

    cache_crosscorr = ~strcmp(cached_dir,'null');  % This is the flag to say if we'll cache results specific to the current spike train.
    if ~exist('psth_option')
        if timevary_PSTH == 0
            psth_option = 0;
        else
            psth_option = 1;
        end
    end

    hashes_of_stims = df_create_stim_cache_file(outputPath,DS);
    hashes_of_spikes = df_create_spike_cache_file(outputPath,DS);
    checksum_CrossCorr = df_checksum('sep',df_load_function_text('df_cal_CrossCorr'),hashes_of_spikes,hashes_of_stims,twindow,smooth_rt,psth_option);
    if isempty(smooth_rt)
        smooth_rt = 41;
    end

    if cache_crosscorr 
        [CSR, CSR_JN, errFlg]= df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_CrossCorr',checksum_CrossCorr,DS,stim_avg,avg_psth,psth,twindow,NBAND);
        save(fullfile(outputPath,'StimResp_crosscorr.mat'), 'CSR');
        save(fullfile(outputPath,'SR_crosscorrJN.mat'), 'CSR_JN');
    else

        [CSR, CSR_JN, errFlg]= df_cal_CrossCorr(DS,stim_avg,avg_psth,...
            psth,twindow,NBAND);
    end
    %     [CSR, CSR_JN, errFlg]= df_cal_CrossCorr(DS,stim_avg,avg_psth,...
    %         psth,twindow,NBAND);

    % Check if df_cal_CrossCorr ends normally
    if errFlg ~= 1        
        % Now call calStrfSep_script to calculate STRF, STRF_JN
        % STRF_JNstd for each tol value
        df_calStrfSep_script;
    end

end

strfFiles = cell(length(Tol_val), 1);
for (k = 1:length(Tol_val))
   fname = [params.outputPath '/strfResult_Tol' num2str(k) '.mat'];
   strfFiles{k} = fname;
end
