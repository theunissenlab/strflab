function spikeTrials = read_spikes_from_file(fileName)

    fid = fopen(fileName, 'r');
    
    spikeStrs = {};
    
    numTrials = 0;
    tl = fgetl(fid);
    lastLineEmpty = 0;
    while (ischar(tl))
       
        [strs] = regexp(tl, ' ', 'split');
        spikeStrs{end+1} = strs;
        lastLineEmpty = isempty(tl);
        numTrials = numTrials + 1;
        tl = fgetl(fid);
    end
    fclose(fid);
    
    if lastLineEmpty
        numTrials = numTrials - 1;
    end
    
    spikeTrials = cell(numTrials, 1);
    for k = 1:numTrials
        
        stimes = [];
        sstr = spikeStrs{k};
        for m = 1:length(sstr)
           
            st = sstr{m};
            if ~isempty(st)
                stimes = [stimes str2double(st)];
            end
        end
        
        spikeTrials{k} = stimes;
    end
    