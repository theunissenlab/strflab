%% Preprocess a series of spike trains into PSTHs
%
%   Input:
%
%       allSpikeTrials: a cell array spike trials, one for each stimulus.
%           Each element contains another cell array, with a length equal
%           to the # of trials. That cell array contains a set of vectors,
%           each vector contains a list of spike times.
%
%       params: Structure of preprocessing parameters:
%
%           .units: time units, either 's' or 'ms', defaults to 's'
%       
%           .stimLengths: length of time in seconds that each stimulus
%               lasts, used to construct PSTH of same length as stimulus
%
%           .binSize: bin size used to create PSTHs, defaults to 0.001s
%
%           .split: if set to 1, even and odd trials are separated out to
%               create two separate PSTHs. These can be used for validation
%               via coherence (see compute_coherence_full.m).
%
%   Output:
%
%       wholeResponse: if params.split=0, this a 1xN vector containing the
%           PSTH, concatenated across trials. if params.split=1, this is a
%           2xN vector, each row containing a PSTH constructed from half
%           the spike trials.
%
%       groupIndex: an 1xN vector representing the group index of each
%           time point, i.e. the stimulus/response pair it belongs to.
%   
%       respInfo: information about the preprocessed response:
%           .responseLengths: length in seconds of each PSTH
%           .numTrials: a vector containing the # of trials for each response
%
%       params: same params as passed into function
%
function [wholeResponse, groupIndex, respInfo, params] = preprocSpikeResponses(allSpikeTrials, params)

    %% set default params
    if ~isfield(params, 'units')
        params.units = 's';
    end
    
    if ~isfield(params, 'binSize')
        params.binSize = 0.001;
    end
    
    if ~isfield(params, 'split')
        params.split = 0;
    end
    
    
    %% turn spike trials into PSTHs    
    totalLength = 0;
    numTrials = zeros(length(allSpikeTrials), 1);
    psthStructs = cell(length(allSpikeTrials), 1);
    respLengths = zeros(1, length(allSpikeTrials));
    for k = 1:length(allSpikeTrials)
       
        pstruct = struct;
        spikeTrials = allSpikeTrials{k};
        for j = 1:length(spikeTrials);        
            %make sure spike times are specified in seconds            
            tmult = 1;
            if strcmp(params.units, 'ms')
                tmult = 1e-3;
            end
            spikeTrials{j} = spikeTrials{j}*tmult;
        end
        numTrials(k) = length(spikeTrials);
        
        stimLen = params.stimLengths(k);
        if ~params.split
            
            pstruct.psth = make_psth(spikeTrials, stimLen, params.binSize);
            psthLen = length(pstruct.psth);
            
        else
            %split spike trials into even and odd sets, create PSTH from
            %each of them
            halfSize = floor(length(spikeTrials) / 2);
            spikeTrials1 = cell(halfSize, 1); %odd trials
            spikeTrials2 = cell(halfSize, 1); %even trials
            
            for j = 1:length(spikeTrials)
                jindx = round(floor(j/2)) + mod(j, 2);
                st = spikeTrials{j};
                if mod(j, 2)                                       
                    spikeTrials1{jindx} = st;
                else
                    spikeTrials2{jindx} = st;
                end
            end
            psth1 = make_psth(spikeTrials1, stimLen, params.binSize);
            psth2 = make_psth(spikeTrials2, stimLen, params.binSize);
            psthLen = length(psth1);
            pstruct.psth = zeros(2, psthLen);
            pstruct.psth(1, :) = psth1;
            pstruct.psth(2, :) = psth2;
        end
        
        respLengths(k) = psthLen * params.binSize;
        totalLength = totalLength + psthLen;
        psthStructs{k} = pstruct;
        
    end
    
    
    %% concatenate PSTHs
    respInfo.responseLengths = respLengths;
    respInfo.numTrials = numTrials;
    
    groupIndex = zeros(1, totalLength);
    
    nPsths = 1 + params.split;
    wholeResponse = zeros(nPsths, totalLength);
    cindx = 1;
    for k = 1:length(psthStructs)
       
        pstruct = psthStructs{k};

        plen = size(pstruct.psth, 2);
        tend = cindx + plen - 1;
        trng = cindx:tend;
        
        wholeResponse(:, trng) = pstruct.psth;        
        groupIndex(trng) = k;
        
        cindx = tend + 1;
        
    end
    
    clear psthStructs;
    