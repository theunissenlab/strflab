function nlinfo = estimate_outputnl_from_spline_cv(linResp, psth, numResample, resampleFraction)

    if nargin < 6
        resampleFraction = 0.75; %how much of the data is used for training
    end

    psmooth = 0.5:0.01:1.0;
    
    trainLen = round(resampleFraction*length(linResp));
        
    validErrs = zeros(numResample, length(psmooth));
    
    for k = 1:numResample
       
        %% create a random training and validation set                        
        trainIndex = randsample(length(linResp), trainLen);
        wholeIndex = 1:length(linResp);
        wholeIndex(trainIndex) = 0;
        validIndex = find(wholeIndex > 0);    
                
        %% fit cubic smoothing splines to dataset with diffent parameters
        for m = 1:length(psmooth)
            outputNL = csaps(linResp(trainIndex), psth(trainIndex), psmooth(m));
            predPsth = fnval(outputNL, linResp(validIndex));
            validErrs(k, m) = sum( (psth(validIndex) - predPsth).^2 );
        end        
        
        [berr, bindx] = min(validErrs(k, :));
        fprintf('Sample %d, best psmooth=%0.2f, err=%0.2f\n', k, psmooth(bindx), berr);
    end
    
    %% find best smoothing parameter
    meanErrs = mean(validErrs, 1);
    [bestErr, bestIndex] = min(meanErrs);
    fprintf('Found best psmooth=%0.2f, err=%0.3f\n', psmooth(bestIndex), bestErr);
    
    %% refit all data using best smoothing parameter
    outputNL = csaps(linResp, psth, psmooth(bestIndex));    
    
    nlinfo = struct;
    nlinfo.outputNL = outputNL;
    nlinfo.psmooth = psmooth(bestIndex);
    