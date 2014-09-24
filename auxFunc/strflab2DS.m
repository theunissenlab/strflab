function DS = strflab2DS(allstim, allresp, groupIndex, outputPath, trainingIndex)

    if nargin < 5
        trainingIndex = 1:length(groupIndex);
    end
    
    allstim = allstim(trainingIndex, :);
    allresp = allresp(trainingIndex);
    groupIndex = groupIndex(trainingIndex);

    groups = unique(groupIndex);
    npairs = length(groups);
    
    DS = cell(npairs, 1);    
    for k = 1:length(groups)
        
        gnum = groups(k);
        
        rng = find(groupIndex == gnum);
        stim = allstim(rng, :);
        resp = allresp(rng);
        
        dfds = struct;
        
        %write spectrogram to intermediate file (direct fit requires this)
        stimfile = fullfile(outputPath, sprintf('df_temp_stim_%d.mat', k));
        outSpectrum = stim';
        save(stimfile, 'outSpectrum');
        clear outSpectrum;        
        dfds.stimfiles = stimfile;
        
        %write response to intermediate file
        rfile = fullfile(outputPath, sprintf('df_temp_resp_%d.mat', k));
        rawResp = resp;
        save(rfile, 'rawResp');
        clear rawResp;
        dfds.respfiles = rfile;
        
        dfds.nlen = size(stim, 1);
        dfds.ntrials = 1;
        
        DS{k} = dfds;
    end
    