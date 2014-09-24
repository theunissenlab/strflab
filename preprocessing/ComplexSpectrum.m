function [s, fo, pg] = ComplexSpectrum(input, increment, winLength,samprate,doFFTShift)
%
% ComplexSpectrum complex spectrum
% 	s = ComplexSpectrum(input, increment, winLength, samprate, doFFTShift)
% 	Compute the complex spectrum of an input signal.
%	Each frame is [winLength]-long and
%	starts [increment] samples after previous frame's start.
%	[doFFTShift] (default value true) controls whether the data is shifted so that
%	it's center is at the start and end of the array (so FFT is more cosine phase)
%	Only zero and the positive frequencies are returned.
%
%   The hamming window in the original code was changed to a gaussian window where
%   winLength is 6 times the standard deviation of the gaussian. 
%    Frederic Theunissen 
%      April 2003.
%
%   8/21/2003
%   Junli modified the followings:
%    1. drop parameter "fftLen"
%    2. add parameter "samprate"
%    3. add two more output "fo" for frequency label
%        "to" for time label.
%    4. handle odd/even winLength 
%    5. dont do fftshift
%
%
%    Junli: Modification on 2/26/2004
%        1. Change time_dur of spectrogram to (length(input)-winLength)/increment +1
%           Before we set time_dur is equal to length(input)/increment +1. 
%                Why? 
%           Because we want to match time_dur of spike data.
%    Junli: Modification on 3/12/2004
%        1. Zero padding at the beginning and end of signal
%

%%% ComplexSpectrum %%%
% Malcolm Slaney's code
%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    doFFTShift = 1;
end

if size(input, 1) > 1,
	input = input';
    
    % Zero padding at the beginning and end of signal
    if rem(winLength, 2)  
       input = [zeros(1,(winLength+1)/2) input zeros(1, (winLength+1)/2-1)];
   else
       input = [zeros(1,winLength/2) input zeros(1, winLength/2)];
   end
end;

inputLength = length(input);

if inputLength < winLength,
	input(winLength) = 0;
	inputLength = winLength;
end;

frameCount = floor((inputLength-winLength)/increment)+1;
%frameCount = floor(inputLength/increment)+1;

%%%%%%%%%%%%%%%%%%%%%
% fftLen - define spectrum's frequency bands
%%%%%%%%%%%%%%%%%%%%%
fftLen = winLength;
%fftLen = max(2^(nextpow2(winLength)+1), fftLen);
%fftLen = max(winLength, fftLen);
%fftLen = winLength;

%%%%%%%%%%%%%%%%%%%%%%%%
% Hamming window 
%a = .54;
%b = -.46;
%wr = sqrt(increment/winLength);
%phi = pi/winLength;
%ws = 2*wr/sqrt(4*a*a+2*b*b)*(a + b*cos(2*pi*(0:winLength-1)/winLength + phi));
%%%%%%%%%%%%%%%%%%%%%%%%

% The hamming window code is replaced by a Gaussian window code - F.
% Theunissen
wx2 = ((1:winLength)-winLength/2).^2;
wvar = (winLength/6)^2;

ws = exp(-0.5*(wx2./wvar));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output "s" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(fftLen, 2)   
    % winLength is odd
    s = zeros((fftLen+1)/2+1, frameCount);
else
    % winLength is even 
    s = zeros(fftLen/2+1, frameCount);
end

pg = zeros(1, frameCount);
for i=1:frameCount
        start = (i-1)*increment + 1;
        last = start + winLength - 1;
        f = zeros(fftLen, 1);
        f(1:winLength) = ws.*input(start:last);
	    pg(i) = std(f(1:winLength));
        
if 0
        if doFFTShift
                f = [f(winLength/2+1:winLength) ; ...
                                     zeros(fftLen-winLength, 1) ; ...
                     f(1:winLength/2)];
        end
end

        specslice = fft(f);
        if rem(fftLen, 2)   % winLength is odd
            
            s(:,i) = specslice(1:((fftLen+1)/2+1));
        else
            
            s(:,i) = specslice(1:(fftLen/2+1));
        end
        %s(:,i) = specslice(1:(fftLen/2+1));
end

% Assign frequency_label
if rem(fftLen, 2)   % winLength is odd
    select = [1:(fftLen+1)/2];
else
    select = [1:fftLen/2+1];
end

fo = (select-1)'*samprate/fftLen;

% assign time_label
to = (1:size(s,2)-1)';
