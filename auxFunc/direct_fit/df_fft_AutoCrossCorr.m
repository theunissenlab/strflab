function [fstim, fstim_spike, stim_spike_JNf] = df_fft_AutoCrossCorr(...
          stim, stim_spike, CSR_JN,TimeLag, NBAND,nstd_val )
%
%  df_fft_AutoCrossCorr(auto,cross,cross_jn,twindow,nband,nstd_val)
%       -- Take FFT of Cross Correlation matrix and smooth them
%       -- Take FFT of JN_Cross Correlation and smooth them
%   Argument: 
%       autofile  - autocorrelation 
%       crossfile - cross correlation
%       jnfile    - JN_cross correlation
%       twindow   - time frame to calcualte auto correlation
%       nband     - size of spatio-domain
%       nstd_val  - value for calculating contour 
%  Return:
%       -- fftautoCorr, fftcrossCorr, and fftJN_CrossCorr  
%             STRFPAK: STRF Estimation Software
% Copyright �2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

% Orignally Written by FET, 2001. 
% Modified by JXZ, 2002.

% ===========================================================
% FFT and smoothing stim_spike matrix 
% ===========================================================
ncorr = size(stim, 1);
nb = NBAND;
nt = 2*TimeLag +1;
nJN = length(CSR_JN);

asize(1) = nt;
asize(2) = nb;
w = hanning(nt);

% Why here?
stim_spike = fliplr(stim_spike);
for ib=1:nb
   stim_spike(ib,:)=stim_spike(ib,:).*w';
end

% ==========================================================
%  Read the JN cross-correlation
% ==========================================================
stim_spike_JN = zeros(nb, nt, nJN);
%w = hanning(nt);
for iJN=1:nJN
    CSR = CSR_JN{iJN};
    for ib=1:nb
       % JUNLI: Oct 21: stim_spike_JN(ib,:,iJN) = df_wrev(CSR(ib,:)) .* w';
       stim_spike_JN(ib,:,iJN) = df_wrev(CSR(ib,:)).*w';
    end
end

clear CSR_JN CSR
% =========================================================
%  The rest part is to calculate fft for all JN values.
%     -- taken from load_stim_spike.m (Written by FET)
% =========================================================
stim_spike_JNf = fft(stim_spike_JN,[],2);
stim_spike_JNmf = mean(stim_spike_JNf,3);
stim_spikef = fft(stim_spike,[],2);

JNv = (nJN-1)*(nJN-1)/nJN;
j = sqrt(-1);
hcuttoff=0;
nf = (nt-1)/2 +1;
for ib=1:nb
	itstart = 1;
	itend = nf;
	below = 0;
	for it=1:nf
		stim_spike_JNvf(ib,it) = JNv*cov(permute(real(stim_spike_JNf(ib,it,:)),[3 2 1]))+ j*JNv*cov(permute(imag(stim_spike_JNf(ib,it,:)),[3 2 1]));
%		rmean = real(stim_spike_JNmf(ib,it));
		rmean = real(stim_spikef(ib,it));
      rstd = sqrt(real(stim_spike_JNvf(ib,it)));
%		imean = imag(stim_spike_JNmf(ib,it));
      imean = imag(stim_spikef(ib,it));
		istd = sqrt(imag(stim_spike_JNvf(ib,it)));
		if abs(rmean) < nstd_val*rstd & abs(imean) < nstd_val*istd
			if itstart == 1 
				itstart = it;
			end
			below = below + 1;
			if below == 3
				itend = it;
			end
		end
		stim_spike_sf(ib,it)= rmean + j*imean;
	end
	for it=1:nf
		if it > itstart 
            expval = exp(-0.5*(it-itstart)^2/(itend-itstart)^2);
            stim_spike_sf(ib,it) = stim_spike_sf(ib,it)*expval;
            stim_spike_JNf(ib,it,:) = stim_spike_JNf(ib,it,:).*expval;
% old way
			%if it < itend
				%stim_spike_sf(ib,it) = stim_spike_sf(ib,it)*(itend-it)/(itend-itstart);
				%stim_spike_JNf(ib,it,:) = stim_spike_JNf(ib,it,:).*((itend-it)/(itend-itstart));
            %else
				%stim_spike_sf(ib,it) = 0.0;
				%for iJN=1:nJN
					%stim_spike_JNf(ib,it,iJN) = 0.0;
                %end
            %end
		end
		if it ~= 1
			stim_spike_sf(ib,nt+2-it) = conj(stim_spike_sf(ib,it));
			stim_spike_JNf(ib,nt+2-it,:) = conj(stim_spike_JNf(ib,it,:));
		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if 0

JNv = (nJN-1)*(nJN-1)/nJN;
j = sqrt(-1);
hcuttoff=0;
nf = (nt-1)/2 +1;
for ib=1:nb
    itstart = 1;
    itend = nf;
    below = 0;

    for it=1:nf
	stim_spike_JNvf(ib,it) = JNv*cov(permute(real(stim_spike_JNf(ib,it,:)),[3 2 1]))+...
         j*JNv*cov(permute(imag(stim_spike_JNf(ib,it,:)),[3 2 1]));
	rmean = real(stim_spikef(ib,it));
	rstd = sqrt(real(stim_spike_JNvf(ib,it)));
	imean = imag(stim_spikef(ib,it));
	istd = sqrt(imag(stim_spike_JNvf(ib,it)));
	if abs(rmean) < nstd_val*rstd & abs(imean) < nstd_val*istd
	    if itstart == 1
		itstart = it;
	    end
	    below = below + 1;
	    if below == 3
		itend = it;
	    end
        end
        stim_spike_sf(ib,it)= rmean + j*imean;
    end

    for it=1:nf
        if it > itstart
            expval = exp(-0.5*(it-itstart)^2/(itend-itstart)^2);
            stim_spike_sf(ib,it) = stim_spike_sf(ib,it)*expval;
            stim_spike_JNf(ib,it,:) = stim_spike_JNf(ib,it,:).*expval;
        end
        if it ~= 1
            stim_spike_sf(ib,nt+2-it) = conj(stim_spike_sf(ib,it));
            stim_spike_JNf(ib,nt+2-it,:) = conj(stim_spike_JNf(ib,it,:));
        end
    end
end

end

% stim_spike in the Fourier Domain is calculated in plot_stim_spike2. 
% This is the smoothed version. Here stim_spike_sf is calculated in
% load_stim_spikejxz.
fstim_spike = stim_spike_sf;

% Reverse fft
if 0
stim_spike_JNm = real(ifft(stim_spike_JNmf,[],2));
stim_spike_s = real(ifft(stim_spike_sf,[],2));
stim_spike_JNs = real(ifft(stim_spike_JNf,[],2));
end

% Choose a window for smoothness
%w = hanning(nt);

% Stim fft
nt2=(nt-1)/2;
fstim=zeros(size(stim));
for i=1:ncorr
    sh_stim = zeros(1,nt);
    w_stim = zeros(1,nt);
    w_stim =  stim(i,:).*w';
    sh_stim(1:nt2+1)=w_stim(nt2+1:nt);
    sh_stim(nt2+2:nt)=w_stim(1:nt2);
    fstim(i,:) = fft(sh_stim);
end


% ===========================================
% END of df_fft_AutoCrossCorr.m
% ===========================================

