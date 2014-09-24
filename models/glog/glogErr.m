function [modelParams, err] = glogErr(modelParams, datIdx)

    global globDat;
    
    [modelParams, predResp] = glogFwd(modelParams, datIdx);
    
    errsq = (rv(globDat.resp(datIdx)) - rv(predResp)).^2;
    err = sum(errsq);
    
    