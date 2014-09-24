function distInfo = spiketrig_dists(psth, predPsth, numTrials, rng)

    %% for spike-triggered distributions
    npsth = psth*numTrials;
    respSpike = [];
    indx = find(npsth > 0);
    for k = indx
        pval = round(npsth(k));
        respSpike = [respSpike repmat(predPsth(k), 1, pval)];
    end
    
    respNoSpike = predPsth(npsth == 0);
    
    %% use density estimations to fit p(x|spike) and p(x|nospike)
    [pxspike, xi, bwpx] = ksdensity(respSpike, rng);
    [pxnospike, xi, bwpxnospike] = ksdensity(respNoSpike, rng);
    
    distInfo.pxspike = pxspike;
    distInfo.pxnospike = pxnospike;
    distInfo.klDist = kl_dist(pxspike, pxnospike, rng);
    