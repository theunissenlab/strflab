function [out,savedStimFft] = stim_error_corr(stim,the_error,maxlag,savedStimFft)
%Returns the error-triggered stimulus; useful for computing gradients.

N = size(stim,1);  %The number of time bins of the stimulus.

if length(the_error) ~=N
    error(['Error doing the error-triggered stimulus correlation: the stim is of size ' num2str(size(stim)) ' and I expected it to have  ' num2str(N) ' rows.']);
end

if ~exist('maxlag','var')
    maxlag = N-1;
end

ERROR = conj(fft(the_error,2^nextpow2(2*N-1))); % do this only once, to save a fft in the loop.
out = zeros(2*maxlag+1,size(stim,2));
if nargout > 1
    if ~exist('savedStimFft','var')
        savedStimFft = fft(stim,2^nextpow2(2*N-1));
        %savedStimFft = savedStimFft';
    end



    %bigc = zeros(2^nextpow2(2*N-1),size(stim,2));
    for jj = 1:size(stim,2)
        STIM = (savedStimFft(:,jj));
        c = real(ifft(ERROR.*STIM));
        %bigc(:,jj) = c;
        out(:,jj) = [c(end-maxlag+1:end); c(1:maxlag+1)];
    end
else
    for jj = 1:size(stim,2)
        STIM = fft(stim(:,jj),2^nextpow2(2*N-1));
        c = real(ifft(ERROR.*STIM));
        out(:,jj) = [c(end-maxlag+1:end); c(1:maxlag+1)];
    end
end
        %size(bigc)
        %out = [bigc(end-maxlag+1:end,:) ; bigc(1:maxlag+1,:)];
        %
        % out = zeros(2*maxlag+1,size(stim,2));
        %  for jj = 1:size(stim,2)
        %      out(:,jj) = xcorr(stim(:,jj),the_error,maxlag);
        %  end

        %{
bigc = real(ifft((ERROR*ones(1,size(stim,2))).*savedStimFft));
out = [bigc(end-maxlag+1:end,:) ; bigc(1:maxlag+1,:)];
%size(ERROR)
%size(savedStimFft)
        %}