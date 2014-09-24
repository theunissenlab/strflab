function psth = make_psth(spikeTrials, stimdur, binsize)

    nbins = round(stimdur/binsize);
    psth = zeros(1, nbins);

    ntrials = length(spikeTrials);

    for k = 1:ntrials

        stimes = spikeTrials{k};
        indx = ((stimes > 0) & (stimes < stimdur));

        stimes = stimes(indx);
        sindxs = round(stimes/binsize) + 1;
        sindxs(sindxs < 1) = 1;
        sindxs(sindxs > nbins) = nbins;
        
        psth(sindxs) = psth(sindxs) + 1;    
    end

    psth = psth / ntrials;
