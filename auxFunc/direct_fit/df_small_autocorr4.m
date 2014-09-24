function CS = df_small_autocorr4(stimval,nband,size_CS,twin,use_fourier)
%  Same as small_autocorr, but without the ntrials multiplication, which
%  can be done outside.
xb = 1;
nband = size(stimval,1);
CS =  zeros(nband/2*(nband+1),2*twin+1);% zeros(size_CS);
N = size(stimval,2);
S = nband;
    time_domain_comp_time = 5e-12 * N * (2*twin +1) * S^2.9; %Time in seconds the time-domain autocorr would take
    fourier_domain_comp_time = 1.15e-8 * S^2 * N * log(N+1); %Time in seconds the fourier-domain autocorr would take
if ~exist('use_fourier','var')
    %use_fourier = 1;
    use_fourier = time_domain_comp_time > fourier_domain_comp_time;
end
if use_fourier

    savedStimFft = fft(stimval',2^nextpow2(2*N-1))';

    for ib1 = 1:nband


        for ib2 = ib1:nband
            c = real(ifft(conj(savedStimFft(ib2,:)).*(savedStimFft(ib1, :))));


            % NEW version of algorithm by using xcorr
            %size(c(end-twin+1:end))
            %size(c(1:twin+1))
            CS(xb, :) = [c(end-twin+1:end) c(1:twin+1)];

            xb = xb +1;
        end            % end of ib2
    end                 % end of ib1
else
    st = stimval';
    for tid = -twin:twin
        onevect = (max(1,tid+1)):(min(N,N+tid));
        othervect = onevect - tid;
        temp = stimval(:,onevect)*st(othervect,:);
        CS(:,twin + tid + 1) = df_uppertrivect(temp);
    end
end
