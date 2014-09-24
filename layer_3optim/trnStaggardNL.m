%% Use a staggard optimization procedure to optimize a linear-nonlinear model
% Input:
%
%   modelParams: a combined linear/nonlinear model (lnl)
%
%   trainingIndex: index of training samples
%
%   metaOptions: meta optimization options
%       .linOptOptions: an options structure applied to the linear
%          optimization part, such as the one returned from trnGradDesc(),
%          trnLARS()
%
%       .nlOptOptions: an options structure applied to the nonlinear
%          optimization, such as the one returned from trnSCG().
%
%       .maxIter: maximum # of iterations
%
%       .infoWindowSize: size in seconds of window used to compute multi-
%           tapered coherence (when finding best iteration)
%
%       .infoFreqCutoff: cutoff in Hz for multi-tapered coherence
%
%       .sampleRate: sample rate in Hz of stimulus/response
%
%       .jackknife: do jacknifing when fitting the output nonlinearity
%
%   validationIndex: index of validation samples
%
% Output:
%
%   modelParamsTrained: trained linear and nonlinear model parameters
%
%
function [modelParamsTrained, metaOptions] = trnStaggardNL(modelParams, trainingIndex, metaOptions, validationIndex)

    global globDat
    
    %% set default parameters
    optDeflt.funcName = 'trnStaggardNL';
    optRange.funcName = {'trnStaggardNL'};

    optDeflt.display = 0;
    optDeflt.maxIter = 100;
    optDeflt.infoWindowSize = 0.250;
    optDeflt.infoFreqCutoff = 90;
    optDeflt.jackknife = 0;
    
    if nargin < 1
        modelParamsTrained = optDeflt;
        return;
    end
        
    %% error checking
    if nargin < 3
        error('You must specify linear/nonlinear model params and optimization options for staggard optimization.\n');    
    end
    
    if ~isfield(metaOptions, 'sampleRate')
        error('You must specify a sample rate for trnStaggardNL!');
    end
    
    if ~strcmp('lnl', modelParams.type)
        error('You must specify a linear-nonlinear model (lnl) to use staggard optimization.\n');
    end
    
    invFuncName = sprintf('%sInvert', modelParams.nlModel.type);
    invPath = which(invFuncName);
    if isempty(invPath)
        error('The nonlinear model must have an inversion function with the name of %s.', invPath);
    end
    
    %% set things up
    origStim = globDat.stim;
    origResp = globDat.resp;
    groupIndex = globDat.groupIdx;
    groups = unique(groupIndex);
        
    linModel = modelParams.linModel;
    nlModel = modelParams.nlModel;
    
    linOpts = metaOptions.linOptOptions;
    nlOpts = metaOptions.nlOptOptions;
    
    linErrs = [];
    nlErrs = [];
    linModels = {};
    nlModels = {};
    
    iter = 0;
    converged = 0;
    
    %% iterate, alternating between linear and nonlinear optimization until done
    while ~converged
        
        %% optimize linear part
        [linModel, linOpts] = strfOpt(linModel, trainingIndex, linOpts);
        
        %produce linear outputs with which to optimize over
        [linModel, linResp] = strfFwd(linModel, trainingIndex);
        linResp(isnan(linResp)) = 0;
        
        %compute error
        [linModel, linErr] = strfErr(linModel, trainingIndex);
        linErrs(end+1) = linErr;
        linModels{end+1} = linModel;
        
        %% optimize nonlinear part
        %set global data to the linear response and PSTH
        strfData(linResp, origResp(trainingIndex), groupIndex(trainingIndex));
        %[nlModel, nlOpts] = strfOpt(nlModel, trainingIndex, nlOpts);
                
        %%[gb, gc] = fit_gompertz_fmin(linResp, origResp(trainingIndex), 1, -abs(randn()), -abs(randn()));        
        
        if metaOptions.jackknife
            fprintf('Using jackknifing to fit output nonlinearity...\n');
            
            [bAll, cAll] = fit_gompertz_fmin(linResp, origResp(trainingIndex), 1, -abs(randn()), -abs(randn()));            
            gindx = groupIndex(trainingIndex);
            [bJN, cJN] = fit_gompertz_fmin_jacknife(linResp, origResp(trainingIndex), gindx, 1, -abs(randn()), -abs(randn()));
            fprintf('   Iteration %d: b=%0.3f +/- %0.4f | c=%0.3f +/- %0.4f\n', iter, mean(bJN), std(bJN), mean(cJN), std(cJN));
            N = length(groups);
            b_pvals = N*bAll - (N-1)*bJN;
            c_pvals = N*cAll - (N-1)*cJN;
            
            bCorrected = mean(b_pvals);
            cCorrected = mean(c_pvals);
            
            %bAvg= mean(bJN);
            %cAvg = mean(cJN);            
            %fprintf('bias corrected: b=%f, c=%f | means: b=%f, c=%f\n', bCorrected, cCorrected, bAvg, cAvg);
            
            bToUse = bCorrected;
            cToUse = cCorrected;
            
        else
            [bToUse, cToUse] = fit_gompertz_fmin(linResp, origResp(trainingIndex), 1, -abs(randn()), -abs(randn()));
        end
        
        nlModel.a = 1;
        nlModel.b = bToUse;
        nlModel.c = cToUse;
        nlModels{end+1} = nlModel;
        
        %compute error
        [nlModel, nlErr] = strfErr(nlModel, 1:length(linResp));
        nlErrs(end+1) = nlErr;
                
        %% invert data for next training iteration
        invFuncCall = sprintf('%sInvert(nlModel, origResp)', nlModel.type);
        [nlModel, invResp] = eval(invFuncCall);
        strfData(origStim, invResp, groupIndex);
        
        iter = iter + 1;
        
        %% evaluate convergence criteria
        if iter >= metaOptions.maxIter
            fprintf('Maximum # of iterations reached.\n');
            converged = 1;
        end
        
    end
    
    %% set global data to original stim/response
    strfData(origStim, origResp, groupIndex);
    
    
    %% find best iteration
    bestScore = -1;
    bestIter = -1; 
    
    for k = 1:metaOptions.maxIter
        
        linModel = linModels{k};
        nlModel = nlModels{k};
        
        lnlModel = lnlInit(linModel, nlModel);
        
        trainRespReal = origResp(trainingIndex);
        [lnlModel, trainRespModel] = strfFwd(lnlModel, trainingIndex);
        trainRespModel(isnan(trainRespModel)) = 0;
                
        validRespReal = origResp(validationIndex);
        [lnlModel, validRespModel] = strfFwd(lnlModel, validationIndex);
        validRespModel(isnan(validRespModel)) = 0;
        
        tiData = compute_coherence_mean(trainRespModel, trainRespReal, metaOptions.sampleRate, metaOptions.infoFreqCutoff, metaOptions.infoWindowSize);
        viData = compute_coherence_mean(validRespModel, validRespReal, metaOptions.sampleRate, metaOptions.infoFreqCutoff, metaOptions.infoWindowSize);
        
        %score = 0.35*tiData.info + 0.65*viData.info;
        score = viData.info;
        if score > bestScore
            bestScore = score;
            bestIter = k;
        end
        
        fprintf('Iter: %d | infoTrain=%f | infoValid=%f | score=%f\n', k, tiData.info, viData.info, score);        
    end
    
    fprintf('Best Model found at iteration %d\n', bestIter);
    
    
    metaOptions.diagnostics.linErrs = linErrs;
    metaOptions.diagnostics.nlErrs = nlErrs;
    metaOptions.diagnostics.linModels = linModels;
    metaOptions.diagnostics.nlModels = nlModels;
    metaOptions.diagnostics.bestIteration = bestIter;
    
    modelParamsTrained = modelParams;
    modelParamsTrained.linModel = linModels{bestIter};
    modelParamsTrained.nlModel = nlModels{bestIter};
    