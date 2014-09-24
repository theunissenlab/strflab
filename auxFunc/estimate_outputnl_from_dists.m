function nlinfo = estimate_outputnl_from_dists(linResp, psth, numTrials, groupIndex)

    numDensityPoints = [150 250 350 450];
    psmooth = [0.90 1.00];
    
    minThreshold = [1e-12];
    
    %% break data into 10 partitions
    groups = unique(groupIndex);
    partitions = cv_partition(10, 1:length(groups));    
    
    nparams = length(numDensityPoints)*length(psmooth)*length(minThreshold);
    perfInfo = cell(nparams, 1);
    
    %% estimate spike-triggered densities    
    cnt = 1;
    for k = 1:length(numDensityPoints)
        for m = 1:length(psmooth)            
            for j = 1:length(minThreshold)
                
                %% use cross-validation to judge goodness of parameters
                npts = numDensityPoints(k);
                ps = psmooth(m);
                mthresh = minThreshold(j);
                
                perrs = zeros(length(partitions), 1);
                kldists = zeros(length(partitions), 1);
                for p = 1:length(partitions)
                
                    %% separate into training and validation sets
                    validationGroups = partitions{p}.validation;
                    trainingGroups = partitions{p}.training;                    
                    validationIndex = findIdx(validationGroups, groupIndex);
                    trainingIndex = findIdx(trainingGroups, groupIndex);
                    
                    trainResp = linResp(trainingIndex);                    
                    trainPsth = psth(trainingIndex);
                    validResp = linResp(validationIndex);
                    validPsth = psth(validationIndex);                    
                                        
                    
                    %% fit distributions on training set
                    distinfo = estimate_spike_dists(trainResp, trainPsth, numTrials, npts, ps, mthresh);
                    outputNL = csaps(distinfo.x, distinfo.pspikex_smooth, 1.0);

                    %% evaluate performance on validation set
                    predPsth = fnval(outputNL, validResp);
                    
                    perrs(p) = sum( (validPsth - predPsth).^2 );
                    
                    %compute KL distance between spike and no-spike distributions
                    predPsthAll = fnval(outputNL, linResp);
                    distInfo = spiketrig_dists(psth, predPsthAll, numTrials, distinfo.x);
                    kldists(p) = distInfo.klDist;
                    
                    %fprintf('\tp=%d, err=%f, kldist=%f\n', p, perrs(p), kldists(p));
                    
                end

                errMean = mean(perrs);
                errDev = std(perrs);
                klMean = mean(kldists);
                klDev = std(kldists);
                
                pinfo = struct;
                pinfo.errMean = errMean;
                pinfo.errDev = errDev;
                pinfo.klMean = klMean;
                pinfo.klDev = klDev;
                pinfo.psmooth = ps;
                pinfo.numDensityPoints = npts;
                pinfo.minThreshold = mthresh;
                
                perfInfo{cnt} = pinfo;
                cnt = cnt + 1;
                
                fprintf('[estimate_outputnl_from_dists] npts=%d, psmooth=%0.2f: err=%0.4f +/- %0.4f, kl=%0.4f +/- %0.4f\n', ...
                        npts, ps, errMean, errDev, klMean, klDev);
                   
            end
        end
    end

    %% find best info
    bestErr = Inf;
    bestKL = 0;
    bestParamsErr = -1;
    bestParamsKL = -1;
    for k = 1:length(perfInfo)       
        pinfo = perfInfo{k};
        if pinfo.errMean < bestErr
            bestErr = pinfo.errMean;
            bestParamsErr = pinfo;
        end
        if pinfo.klMean > bestKL
            bestKL = pinfo.klMean;
            bestParamsKL = pinfo;
        end
    end
    
    
    %% debugging
    %{
    qvals = zeros(length(perfInfo), 2);
    for k = 1:length(perfInfo)
        qvals(k, 1) = perfInfo{k}.errMean;
        qvals(k, 2) = perfInfo{k}.klMean;
    end 
    
    figure; hold on;
    plot(qvals(:, 1), qvals(:, 2), 'ko');
    axis tight; xlabel('msqerr'); ylabel('kl dist');
    %}
    
    %% fit on all data using best parameters
    distinfoErr = estimate_spike_dists(linResp, psth, numTrials, ...
                                       bestParamsErr.numDensityPoints, bestParamsErr.psmooth, bestParamsErr.minThreshold);
    outputNLErr = csaps(distinfoErr.x, distinfoErr.pspikex_smooth, 1.0);
    
    distinfoKL = estimate_spike_dists(linResp, psth, numTrials, ...
                                       bestParamsKL.numDensityPoints, bestParamsKL.psmooth, bestParamsKL.minThreshold);
    outputNLKL = csaps(distinfoKL.x, distinfoKL.pspikex_smooth, 1.0);
    
    fprintf('Found lowest error with npts=%d, psmooth=%0.2f: err=%0.4f +/- %0.4f, kl=%0.1f +/- %0.1f\n', ...
            bestParamsErr.numDensityPoints, bestParamsErr.psmooth, bestParamsErr.errMean, bestParamsErr.errDev, bestParamsErr.klMean, bestParamsErr.klDev);
    fprintf('Found largest KL distance with npts=%d, psmooth=%0.2f: dist=%0.4f +/- %0.4f, kl=%0.1f +/- %0.1f\n', ...
            bestParamsKL.numDensityPoints, bestParamsKL.psmooth, bestParamsKL.errMean, bestParamsKL.errDev, bestParamsKL.klMean, bestParamsKL.klDev);
        
    nlinfo = struct;
    nlinfo.err = distinfoErr;
    nlinfo.err.outputNL = outputNLErr;
    
    nlinfo.kl = distinfoKL;
    nlinfo.kl.outputNL = outputNLKL;
    
