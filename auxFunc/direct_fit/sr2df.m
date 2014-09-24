%% Convert stim/response data structure (from preprocess.sound.m) to the DS
%% structure that direct_fit uses
%
%   Input:
%       srData: a stim/response data structure produced by
%         preprocess_sound.m
%
%       outputPath: an arbitrary directory where intermediate files are
%         written, used to load stim/responses by direct_fit
%
%   Output:
%
%       DS: a direct_fit friendly dataset representation
%
%   Author: Mike Schachter (mike.schachter@gmail.com)
%
function DS = sr2df(srData, outputPath)
    
    srDatasets = srData.datasets;
    
    npairs = length(srDatasets);
    DS = cell(npairs, 1);    
    for k = 1:npairs
        
        srds = srDatasets{k};
        
        dfds = struct;
        
        %write spectrogram to intermediate file (direct fit requires this)
        stimfile = fullfile(outputPath, sprintf('df_temp_stim_%d.mat', k));
        outSpectrum = srds.stim.tfrep.spec;
        save(stimfile, 'outSpectrum');
        clear outSpectrum;        
        dfds.stimfiles = stimfile;
        
        %write response to intermediate file
        rfile = fullfile(outputPath, sprintf('df_temp_resp_%d.mat', k));
        rawResp = srds.resp.rawSpikeIndicies;
        save(rfile, 'rawResp');
        clear rawResp;
        dfds.respfiles = rfile;
        
        dfds.nlen = size(srds.stim.tfrep.spec, 2);
        dfds.ntrials = length(srds.resp.rawSpikeIndicies);
        
        DS{k} = dfds;
    end
    