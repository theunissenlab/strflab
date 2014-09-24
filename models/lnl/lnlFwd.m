function [modelParams, resp] = lnlFwd(modelParams, datIdx)

    global globDat;
    
    origStim = globDat.stim;
    origResp = globDat.resp;
    groupIndex = globDat.groupIdx;
    
    %% compute linear part of response
    [lp, linResp] = strfFwd(modelParams.linModel, datIdx);
    
    %% compute nonlinear part of response
    strfData(linResp, origResp(datIdx), groupIndex(datIdx));    
    [nlp, resp] = strfFwd(modelParams.nlModel, 1:length(linResp));
    
    %% reset global variables
    strfData(origStim, origResp, groupIndex);
    
    