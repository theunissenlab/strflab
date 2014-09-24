function distinfo = estimate_spike_dists_cv(linResp, psth, numTrials, numResample, resampleFraction)

    numDensityPoints = 50:25:300;
    
    validErrs = zeros(length(numDensityPoints), numResample);

    %% find the best number of density points by cross validation
    for k = 1:length(numDensityPoints)
        
        minX = min(linResp);
        maxX = max(linResp);
        xRng = linspace(minX, maxX, numDensityPoints(k));
        
        for m = 1:length(numResample)        
            trainLen = round(resampleFraction*length(linResp));
            trainIndex = randsample(length(linResp), trainLen);
            wholeIndex = 1:length(linResp);
            wholeIndex(trainIndex) = 0;
            validIndex = find(wholeIndex > 0);    

            xSamp = linResp(trainIndex);
            psthSamp = psth(trainIndex);

            pdata = esd_fit_dists(xSamp, psthSamp, xRng, numTrials);

            %% compute cubic spline approximation and error on validation set        
            fitP = pdata.pspikex;
            fitP(isnan(fitP)) = 0;
            outputNL = csaps(xRng, fitP, 1.0);
            predPsth = fnval(outputNL, linResp(validIndex));
            validErrs(k, m) = sum( (cv(psth(validIndex)) - cv(predPsth)).^2 );            
        end
         
        meanErr = mean(validErrs(k, :));
        meanStd = std(validErrs(k, :));
        fprintf('Sample %d, numDensityPoints=%d, avg err=%0.2f +/- %0.3f\n', k, numDensityPoints(k), meanErr, meanStd);
    end
    
    %% get the best # of density points
    meanErrs = mean(validErrs, 2);
    [bestErr, bestIndex] = min(meanErrs);
    fprintf('\nFound best numDensityPoints=%d, err=%0.3f\n', numDensityPoints(bestIndex), bestErr);
    
    %% resample final density using all data with best # of points
    minX = min(linResp);
    maxX = max(linResp);
    xRng = linspace(minX, maxX, numDensityPoints(bestIndex));
    
    pdata = esd_fit_dists(xSamp, psthSamp, xRng, numTrials);
    fitP = pdata.pspikex;
    fitP(isnan(fitP)) = 0;
    outputNL = csaps(xRng, fitP, 1.0);

    %% create output structure
    distinfo.x = xRng;
    distinfo.px = pdata.px;
    distinfo.pxspike = pdata.pxspike;
    distinfo.pxnospike = pdata.pxnospike;
    distinfo.pspikex = pdata.pspikex;
    distinfo.pnospikex = pdata.pnospikex;
        
    distinfo.outputNL = outputNL;
    
end
    

function pdata = esd_fit_dists(x, psth, rng, numTrials)

    %% construct p(x)
    [px, xi, bwpx] = ksdensity(x, rng, 'function', 'pdf');
    
    %% construct p(x|spike)
    npsth = psth*numTrials;
    xSpike = [];
    indx = find(npsth > 0);
    for m = indx
        pval = round(npsth(m));
        xSpike = [xSpike repmat(x(m), 1, pval)];
    end    
    [pxspike, xi, bwpxspike] = ksdensity(xSpike, rng, 'function', 'pdf');        

    %% construct p(x|nospike)
    xNoSpike = x(npsth == 0);
    [pxnospike, xi, bwpxnospike] = ksdensity(xNoSpike, rng, 'function', 'pdf');

    pspike = mean(psth);
    pnospike = 1 - pspike;

    %% construct p(spike|x) and p(nospike|x)
    pspikex = (pxspike*pspike) ./ px;    
    pnospikex = (pxnospike*pnospike) ./ px;
    
    %% construct a histogram of x, artificially extending it to the size of rng
    [nx, cx] = hist([min(rng) max(rng) cv(x)], length(rng));
    
    sampleThresh = 2;
    zindx = nx < sampleThresh;        
    
    px(zindx) = NaN;
    pxspike(zindx) = NaN;
    pxnospike(zindx) = NaN;
    pspikex(zindx) = NaN;
    pnospikex(zindx) = NaN;
    
    pdata = struct;
    pdata.px = px;
    pdata.pxspike = pxspike;
    pdata.pxnospike = pxnospike;
    pdata.pspike = pspike;
    pdata.pnospike = pnospike;
    pdata.pspikex = pspikex;
    pdata.pnospikex = pnospikex;
end

