function nlinfo = estimate_outputnl_from_spline(linResp, psth, numTrials, groupIndex)
    
    psmooth = 0.5:0.05:1.0;
    
    numDensityPoints = 150:50:300;    

    %% break data into 10 partitions
    groups = unique(groupIndex);
    partitions = cv_partition(10, groups);
    
    perfInfo = cell(length(psmooth), 1);
        
    for k = 1:length(psmooth)
        
        %% use cross-validation to judge goodness of parameters
        ps = psmooth(k);
        
        perrs = zeros(length(partitions), 1);
        klDists = zeros(length(partitions), 1);
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
        
            outputNL = csaps(trainResp, trainPsth, ps);

            predPsth = fnval(outputNL, validResp);
            perrs(p) = sum( (validPsth - predPsth).^2 );
            
            %compute KL distance between spike and no-spike distributions
            % (take mean across different # of sampling points)
            kd = zeros(length(numDensityPoints), 1);
            for m = 1:length(kd)
                xinc = (max(linResp) - min(linResp)) / numDensityPoints(m);
                xrng = min(linResp):xinc:max(linResp);
                predPsthAll = fnval(outputNL, linResp);
                distInfo = spiketrig_dists(psth, predPsthAll, numTrials, xrng);
                kd(m) = distInfo.klDist;
            end
            klDists(p) = mean(kd);
            %fprintf('kl=%0.1f +/- %0.1f\n', mean(kd), std(kd));
        end
        
        %compute least square error
        errMean = mean(perrs);
        errStd = std(perrs);
        
        %compute KL distance
        klMean = mean(klDists);
        klStd = std(klDists);
        
        fprintf('[estimate_outputnl_from_spline] psmooth=%0.2f: err=%0.4f +/- %0.4f, kl=%0.1f +/- %0.1f\n', ps, errMean, errStd, klMean, klStd);
        
        pinfo = struct;
        pinfo.errMean = errMean;
        pinfo.errStd = errStd;
        pinfo.klMean = klMean;
        pinfo.klStd = klStd;
        pinfo.psmooth = ps;
        
        perfInfo{k} = pinfo;
    end

    %% find best parameters    
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

    fprintf('Found minimum error with psmooth=%0.2f: err=%f\n', bestParamsErr.psmooth, bestParamsErr.errMean);    
    fprintf('Found largest KL distance with psmooth=%0.2f: dist=%f\n', bestParamsKL.psmooth, bestParamsKL.klMean);    
    
    
    %% fit entire dataset with best params
    outputNLErr = csaps(linResp, psth, bestParamsErr.psmooth);
    outputNLKL = csaps(linResp, psth, bestParamsKL.psmooth);
    
    nlinfo = struct;
    nlinfo.err = struct;
    nlinfo.err.psmooth = bestParamsErr.psmooth;
    nlinfo.err.outputNL = outputNLErr;
    
    nlinfo.kl = struct;
    nlinfo.kl.psmooth = bestParamsKL.psmooth;
    nlinfo.kl.outputNL = outputNLKL;
    