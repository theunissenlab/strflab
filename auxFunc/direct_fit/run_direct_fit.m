function strfs = run_direct_fit(srData, params)

    if nargin < 2
       fprintf('Using default parameters.\n');    
       params = struct;
    end
    
    %% set default parameters
    flds = {'NBAND', 'Tol_val', 'setSep', 'TimeLag', 'TimeLagUnit', ...
            'timevary_PSTH', 'smooth_rt', 'ampsamprate', 'respsamprate', ...
            'outputPath', 'use_alien_space', 'alien_space_file'};
    dtols = [0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
    dvals = {0, dtols, 0, 100, 'msec', 0, 41, 1000, 1000, tempdir(), 0, ''};        
    params = check_fields(params, flds, 'Direct field requires params.%s!\n', dvals);
    
    fprintf('Using output directory %s\n', params.outputPath);
    
    %% convert stim/response representation into direct fit representation    
    DS = sr2df(srData, params.outputPath);
    
    %% set params from stim/response data
    params.NBAND = srData.nStimChannels;
    params.ampsamprate = srData.stimSampleRate;
    params.respsamprate = srData.respSampleRate;
    params.DS = DS;
    
    fprintf('time varying: %d\n', params.timevary_PSTH);
    
    strfFiles = direct_fit(params);
    
    strfs = cell(length(strfFiles), 1);
    for k = 1:length(strfs)       
        svars = load(strfFiles{k});
        strfs{k} = svars.STRF_Cell;
        clear svars;
    end
