function distinfo = estimate_spike_dists(linResp, psth, numTrials, numDensityPoints, splineSmooth, minThreshold)

    if nargin < 4
        numDensityPoints = 250;
    end
    if nargin < 5
        splineSmooth = 0.60;
    end
    if nargin < 6
        minThreshold = 1e-12;
    end
    
    
    %% construct a range of x to evaluate densities    
    minLin = min(linResp);
    maxLin = max(linResp);
    linWidth = (maxLin - minLin);

    linEps = 0.05*linWidth;    
    minLin = minLin - linEps;
    maxLin = maxLin + linEps;
    
    nSampPoints = min(numDensityPoints, length(linResp));
    linInc = linWidth / nSampPoints;
    linRng = minLin:linInc:maxLin;        

    
    %% compute p(spike) and p(nospike)
    pspike = mean(psth);
    pnospike = 1 - pspike;
    
    
    %% construct p(x) and smooth it
    [px, xi, bwpx] = ksdensity(linResp, linRng);     
    [px_cdf, xi, bwpx] = ksdensity(linResp, linRng, 'function', 'cdf');     
    px_smooth = smooth_pdens(linRng, px, splineSmooth);
        
    
    %% construct p(x|spike) and smooth it
    npsth = psth*numTrials;
    linRespSpike = [];
    indx = find(npsth > 0);
    for k = indx
        pval = round(npsth(k));
        linRespSpike = [linRespSpike repmat(linResp(k), 1, pval)];
    end    
    [pxspike, xi, bwpxspike] = ksdensity(linRespSpike, linRng, 'function', 'pdf');        
    [pxspike_cdf, xi, bwpxspike] = ksdensity(linRespSpike, linRng, 'function', 'cdf');        
    %[bwpxspike, pxspike, pxspike_rng, pxspike_cdf] = kde(rv(linRespSpike), nSampPoints, minLin, maxLin);
    pxspike_smooth = smooth_pdens(linRng, pxspike, splineSmooth);
    
    
    %% construct p(x|nospike) and smooth it
    linRespNoSpike = linResp(npsth == 0);
    [pxnospike, xi, bwpxnospike] = ksdensity(linRespNoSpike, linRng);
    [pxnospike_cdf, xi, bwpxnospike] = ksdensity(linRespNoSpike, linRng, 'function', 'cdf');
    %[bwpxnospike, pxnospike, pxnospike_rng, pxnospike_cdf] = kde(rv(linRespNoSpike), nSampPoints, minLin, maxLin);
    pxnospike_smooth = smooth_pdens(linRng, pxnospike, splineSmooth);
    
    
    %% compute p(spike|x)    
    pspikex = (pxspike*pspike) ./ px;
    pspikex_smooth = (pxspike_smooth*pspike) ./ px_smooth;
    pspikex_smooth(pspikex_smooth < minThreshold) = 0;
    pspikex_smooth(isinf(pspikex_smooth)) = 0;
    pspikex_smooth(isnan(pspikex_smooth)) = 0;
    
    
    %% compute p(nospike|x)    
    pnospikex = (pxnospike*pnospike) ./ px;
    pnospikex_smooth = (pxnospike_smooth*pnospike) ./ px_smooth;
    pnospikex_smooth(pnospikex_smooth < minThreshold) = 0;
    pnospikex_smooth(isinf(pnospikex_smooth)) = 0;
    pnospikex_smooth(isnan(pnospikex_smooth)) = 0;

    
    %% set up output structure
    distinfo = struct;
    distinfo.x = linRng;
    
    distinfo.px = px;
    distinfo.px_smooth = px_smooth;    
    distinfo.px_cdf = px_cdf;
    
    distinfo.pxspike = pxspike;
    distinfo.pxspike_smooth = pxspike_smooth;
    distinfo.pxspike_cdf = pxspike_cdf;
    
    distinfo.pxnospike = pxnospike;
    distinfo.pxnospike_smooth = pxnospike_smooth;
    distinfo.pxnospike_cdf = pxnospike_cdf;
    
    distinfo.pspikex = pspikex;
    distinfo.pspikex_smooth = pspikex_smooth;
    
    distinfo.pnospikex = pnospikex;
    distinfo.pnospikex_smooth = pnospikex_smooth;

end


function px_smooth = smooth_pdens(xrng, px, psmooth)
    pxFunc = csaps(xrng, px, psmooth);
    px_smooth = fnval(pxFunc, xrng);      
    px_smooth(px_smooth < 0) = 0;
    px_smooth(isnan(px_smooth)) = 0;
    px_smooth(isinf(px_smooth)) = 0;
end
    
