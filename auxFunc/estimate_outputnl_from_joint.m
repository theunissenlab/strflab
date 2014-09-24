function [nlinfo, params] = estimate_outputnl_from_joint(linResp, psth, numTrials, params)

    if nargin < 4
        params = struct;
    end
    
    %% set default parameters
    if ~isfield(params, 'smoothingWidth')
        params.smoothingWidth = 0.003;
    end
    if ~isfield(params, 'numJointPointsX')
        params.numJointPointsX = 75;
    end    
    if ~isfield(params, 'numJointPointsY')
        params.numJointPointsY = 25;
    end
    if ~isfield(params, 'numDensityPoints')
        params.numDensityPoints = 250;
    end
    if ~isfield(params, 'splineSmooth')
        params.splineSmooth = 0.85;
    end
    
    %% constuct range of points to sample at
    minLin = min(linResp);
    maxLin = max(linResp);
    linWidth = (maxLin - minLin);
    linEps = 0.05*linWidth;    
    minLin = minLin - linEps;
    maxLin = maxLin + linEps;
    
    %nSampPoints = min(2^7, length(linResp));
    nSampPoints = min(params.numDensityPoints, length(linResp));
    linInc = linWidth / nSampPoints;
    linRng = minLin:linInc:maxLin;
    
    
    %% estimate p(x) density and it's CDF
    [px, xi, bwpx] = ksdensity(linResp, linRng);     
    [px_cdf, xi, bwpx] = ksdensity(linResp, linRng, 'function', 'cdf');     
    %[bwpx, px,px_rng, px_cdf] = kde(rv(linResp), nSampPoints, minLin, maxLin);
    
    
    %% construct p(x|spike)
    npsth = psth*numTrials;
    linRespSpike = [];
    indx = find(npsth > 0);
    for k = indx
        pval = round(npsth(k));
        linRespSpike = [linRespSpike repmat(linResp(k), 1, pval)];
    end    
    [pxspike, xi, bwpxspike] = ksdensity(linRespSpike, linRng, 'function', 'pdf');        
    [pxspike_cdf, xi, bwpxspike] = ksdensity(linRespSpike, linRng, 'function', 'cdf');        

    
    %% smooth the psth
    gwidth = params.smoothingWidth;
    grng = -0.075:0.001:0.075;
    gconv = exp(-(grng / gwidth).^2);

    spsth = conv(psth, gconv, 'same');
    spsth = (spsth / max(spsth))*max(psth);
    
    
    %% now produce joint distribution between smoothed PSTH and linear outputs
    data = [rv(linResp) rv(spsth)];
    
    fitMinVal = min(linRng);
    fitMaxVal = max(linRng);
    gridxInc = (fitMaxVal - fitMinVal) / params.numJointPointsX;
    gridx = fitMinVal:gridxInc:fitMaxVal;
    gridyInc = 1 / params.numJointPointsY;
    gridy = 0:gridyInc:1;
    [pxy, bw] = kde2d_nn(data, gridx, gridy);

    
    %% produce predicted output nonlinearty E[psth|linResp] from joint
    psthPreds = zeros(1, size(pxy, 2));
    for k = 1:length(psthPreds)
        pvals = pxy(:, k);        
        
        xval = gridx(k);
        pxMax = xval + gridxInc;
        pxMin = xval - gridxInc;
        
        pxval = mean(px(linRng >= pxMin & linRng <= pxMax));
        
        psthPreds(k) = ((cv(pvals)*rv(gridy)) / pxval);
    end
    
    
    %% identify X range to fit
    lowXIndx = max(find(pxspike_cdf <= 0.0001));
    lowX = linRng(lowXIndx);
    
    %% set predictions outside of this range to 0
    psthPreds(gridx <= lowX) = 0;
  
    
    %% smooth the output NL prediction with a cubic spline
    outputNL = csaps(gridx, psthPreds, params.splineSmooth);
    
    
    %% set up output parameters
    nlinfo = struct;
    
    nlinfo.pxy = pxy;
    nlinfo.px = px;
    nlinfo.px_cdf = px_cdf;
    nlinfo.pxy_data = data;
    nlinfo.rawNL = psthPreds;
    nlinfo.x = gridx;
    nlinfo.y = gridy;
    
    nlinfo.outputNL = outputNL;
    
    