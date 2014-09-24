function [CSR_JN] = df_internal_cal_CrossCorr(stimval,psthval,twin,do_fourier);
nband = size(stimval,1);
CSR_JN = zeros(nband,2*twin+1);
N = length(psthval);
td_time = 2.5e-8*nband*N*(1+2*twin); %time in s to calculate using a time-domain algorithm
fd_time = 2e-7*N*log(N+1)*nband; %time in s to calculate using a Fuorier-domain algorithm
if ~exist('do_fourier','var')
    do_fourier = fd_time < td_time;
end
if do_fourier
    for ib1 = 1:nband
        CSR_JN(ib1,:) = xcorr(stimval(ib1, :), psthval, twin);
    end
else
    pt = psthval';
    for tid = -twin:twin
        onevect = (max(1,tid+1)):(min(N,N+tid));
        othervect = onevect - tid;
        temp = stimval(:,onevect)*pt(othervect);
        CSR_JN(:,tid+twin+1) = temp;
    end
end
