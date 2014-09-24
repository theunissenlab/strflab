function partitions = cv_partition(nDesiredFolds, groups)

    nSampsPerFold = round(length(groups) / nDesiredFolds);
    nFolds = floor(length(groups)/nSampsPerFold) + mod(length(groups), nSampsPerFold);
    partitions = cell(nFolds, 1);
    
    for k = 1:nFolds        
        b = 1:nSampsPerFold;
        g = (k-1)*nSampsPerFold + b;
        g(g > length(groups)) = [];
        
        tg = groups;
        for m = 1:length(g)
            tg(tg == g(m)) = [];
        end
        
        pstruct = struct;
        pstruct.training = tg;
        pstruct.validation = g;
        
        partitions{k} = pstruct;        
    end
    
    