function [allstim, allresp, groupIndex] = srdata2strflab(srData)

    pairCount = length(srData.datasets);
    totalStimLength = 0;
    totalRespLength = 0;
    
    %tally up total lengths
    for k = 1:pairCount
        ds = srData.datasets{k};
        totalStimLength = totalStimLength + size(ds.stim.tfrep.spec, 2);
        totalRespLength = totalRespLength + length(ds.resp.psth);
    end
    
    allstim = zeros(totalStimLength, srData.nStimChannels);
    allresp = zeros(1, totalRespLength);
    groupIndex = zeros(1, totalRespLength);
    
    %set up matricies
    currentIndex = 1;
    for k = 1:pairCount
        ds = srData.datasets{k};
        
        stim = ds.stim.tfrep.spec;
        resp = ds.resp.psth;
        
        if size(stim, 2) ~= length(resp)
            error('Stim and response lengths are not the same for dataset %d!\n', k);
        end
        eIndx = currentIndex + length(resp) - 1;
        rng = currentIndex:eIndx;
        
        allstim(rng, :) = stim';
        allresp(rng) = resp;
        groupIndex(rng) = k;
        
        currentIndex = eIndx+1;
    end
    
    