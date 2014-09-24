function [modelParams, err] = gompErr(modelParams, datIdx)

    global globDat;
    
    [modelParams, predResp] = gompFwd(modelParams, datIdx);
    
    errsq = (rv(globDat.resp(datIdx)) - rv(predResp)).^2;
    err = sum(errsq);
    
    