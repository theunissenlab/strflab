function distinfo = estimate_spike_dists_cv2(linResp, psth, numTrials, numDensityPoints, numResample, resampleFraction)

    if nargin < 6
        resampleFraction = 0.75; %how much of the data is resampled
    end
    
    resampleCount = round(resampleFraction*length(linResp));
    
    minX = min(linResp);
    maxX = max(linResp);
    xRng = linspace(minX, maxX, numDensityPoints);
    
    psmooth = 0.5:0.01:1.0;
    validErrs = zeros(numResample, length(psmooth));
   
    px = zeros(numResample, length(xRng));
    pxspike = zeros(numResample, length(xRng));
    pxnospike = zeros(numResample, length(xRng));
    pspikex = zeros(numResample, length(xRng));
    pnospikex = zeros(numResample, length(xRng));
    
    pspike = zeros(numResample, 1);
    pnospike = zeros(numResample, 1);
    
    fprintf('Samples have %d points\n', resampleCount);
    
    for k = 1:numResample
        
        trainLen = round(resampleFraction*length(linResp));
        trainIndex = randsample(length(linResp), trainLen);
        wholeIndex = 1:length(linResp);
        wholeIndex(trainIndex) = 0;
        validIndex = find(wholeIndex > 0);    
    
        xSamp = linResp(trainIndex);
        psthSamp = psth(trainIndex);
        
        pdata = esd_fit_dists(xSamp, psthSamp, xRng, numTrials);
        
        px(k, :) = pdata.px;
        pxspike(k, :) = pdata.pxspike;
        pxnospike(k, :) = pdata.pxnospike;
        
        pspike(k) = pdata.pspike;
        pnospike(k) = pdata.pnospike;
        pspikex(k, :) = pdata.pspikex;
        pnospikex(k, :) = pdata.pnospikex;
        
        %% compute cubic spline approximation and error on validation set        
        fitP = pspikex(k, :);
        fitP(isnan(fitP)) = 0;
        for m = 1:length(psmooth)
            outputNL = csaps(xRng, fitP, psmooth(m));
            predPsth = fnval(outputNL, linResp(validIndex));
            validErrs(k, m) = sum( (cv(psth(validIndex)) - cv(predPsth)).^2 );
        end
        validErrs(k, :)
        [berr, bindx] = min(validErrs(k, :));
        fprintf('Sample %d, %d points, best psmooth=%0.2f, err=%0.2f\n', k, trainLen, psmooth(bindx), berr);
    end
    
    %% find best smoothing parameter
    meanErrs = mean(validErrs, 1);
    [bestErr, bestIndex] = min(meanErrs);
    fprintf('\nFound best psmooth=%0.2f, err=%0.3f\n', psmooth(bestIndex), bestErr);
    
    %% fit final smoother on entire data
    pdata = esd_fit_dists(linResp, psth, xRng, numTrials);
    fitP = pdata.pspikex;
    fitP(isnan(fitP)) = 0;
    outputNL = csaps(xRng, fitP, psmooth(bestIndex));

    %% create output structure
    distinfo.x = xRng;
    distinfo.px = px;
    distinfo.pxspike = pxspike;
    distinfo.pxnospike = pxnospike;
    distinfo.pspikex = pspikex;
    distinfo.pnospikex = pnospikex;
    
    distinfo.final = pdata;
    
    distinfo.psmooth = psmooth(bestIndex);
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

