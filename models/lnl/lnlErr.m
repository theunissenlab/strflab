function [modelParams, err] = lnlErr(modelParams, datIdx)

    global globDat;
    
    [modelParams, predResp] = lnlFwd(modelParams, datIdx);
    
    errsq = (globDat.resp(datIdx) - predResp).^2;
    err = sum(errsq);
    
    