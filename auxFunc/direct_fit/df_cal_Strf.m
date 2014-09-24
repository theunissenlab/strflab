function [strfH, strfHJN, strfHJN_std] = df_cal_Strf(fstim, fstim_spike,...
    stim_spike_JNf, stim_size, stim_spike_size, stim_spike_JNsize,...
    nb, nt, nJN, tol, save_flag)
%
%  [strfH, strfHJN, strfHJN_std] = df_cal_Strf(fstim, fstim_spike,...
%     stim_spike_JNf, stim_size, stim_spike_size, stim_spike_JNsize,...
%     nb, nt, nJN, tol)
%      -- Calculate strf, JackKnifed strf and contour JN strf
%     Input:
%         fstim: FFT of stim auto-correlation
%         fstim_spike: FFT of stim_spike cross-correlation
%         stim_spike_JNf: FFT of JackKnifed stim_spike cross-correlation
%         stim_size: size of stim auto-correlation
%         stim_spike_size: size of stim_spike cross-correlation
%         stim_spike_JNsize: size of JN stim_spike cross-correlation
%         nb: length of spatio domain
%         nt: length of time domain
%         nJN: num of JackKnife case
%         tol: one tol. value
%         save_flag: the flag to save the intermediate result for specific tol
%                    The default value is 0 (no save), 1 otherwise.
%     Output:
%         strfH: the estimated strf
%         strfHJN: the estimated JackKnifed version strf
%         strfHJN_std: the estimated JN strf
%
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

% Modified by JXZ
% 7/12/2005
%     1. Add 'durs' -- stimuli's duration
%     2. Recalculate 'cross_vectJN' for 'forwardJN_std'
%     3. Recalculate 'forwardJN_std'
%

global DF_PARAMS

% =======================================================
% Check if we have all input
% =======================================================
if ~exist('save_flag')
    save_flag = 0;
end

% =======================================================
%  JXZ: 7/12/2005
%  Try to get stimuli's durations from variable 'DS'
% =======================================================
DS = DF_PARAMS.DS;
for jj = 1:length(DS)
    durs(jj) = DS{jj}.nlen;
end

% ========================================================
% Forward Filter - The algorithm is from FET's filters2.m
% ========================================================
nf = (nt-1)/2 +1;

% Allocate space for all arrays
clear j;
stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
cross_vect = zeros(nb,1);
cross_vectJN = zeros(nJN,nb);
h = zeros(1,nb);
hJN = zeros(nJN, nb);
ffor=zeros(stim_spike_size);
fforJN=zeros(stim_spike_JNsize);
strfH =zeros(stim_spike_size);
strfHJN=zeros(stim_spike_JNsize);
cums=zeros(nf,nb+1);
ranktest=zeros(1,nf);
stimnorm=zeros(1,nf);
clear j;

test=ones(nb,nb);

% ========================================================
% Find the maximum norm of all the matrices
% ========================================================
for iff=1:nf
    % I great big thanks to Georg for the next lines, which improve speed:
    stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
    stim_mat(find(tril(test)))=conj(fstim(:,iff));
    stim_mat=stim_mat-diag(diag(stim_mat))+(stim_mat');
    %
    %         nc = 1;
    %         for fb_indx=1:nb
    %             for fb_indx2=fb_indx:nb
    %                 stim_mat(fb_indx,fb_indx2) = fstim(nc,iff);
    %                 if fb_indx ~= fb_indx2
    %                     stim_mat(fb_indx2,fb_indx) = conj(fstim(nc,iff));
    %                 end
    %                 nc = nc +1;
    %             end
    %         end
    stimnorm(1,iff)=norm(stim_mat);
end

% One way of doing it
ranktol=tol*max(stimnorm);

% Do the matrix inversion for each frequency
for iff=1:nf
    % Stuff stim matrix and cross-correlation vectors
    stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
    stim_mat(find(tril(test)))=conj(fstim(:,iff));
    stim_mat=stim_mat-diag(diag(stim_mat))+(stim_mat');

    nc = 1;
    for fb_indx=1:nb
        %         for fb_indx2=fb_indx:nb
        %             stim_mat(fb_indx,fb_indx2) = fstim(nc,iff);
        %             if fb_indx ~= fb_indx2
        %                 stim_mat(fb_indx2,fb_indx) = conj(fstim(nc,iff));
        %             end
        %             nc = nc +1;
        %         end
        % cross_vect(i) = conj(fstim_spike(i,iff));
        cross_vect(fb_indx) = fstim_spike(fb_indx,iff);
        for iJN=1:nJN
            % =======================================================
            %  JXZ: 7/12/2005
            % cross_vectJN(iJN,i) = conj(stim_spike_JNf(i,iff,iJN));
            % =======================================================
            cross_vectJN(iJN,fb_indx) = stim_spike_JNf(fb_indx,iff,iJN);
            %Normalized by the duration of each individual stim
            % This code does not make sense to me becuase the stim_spike_JN
            % is already the delete 1 (correctly normalized).
%             nstim = length(durs);
%             jn2_index = 1+mod((1:(nstim-1))+iJN,nstim);  %All but the stim after the current stim; used in double-jackknife
%             n_jn2 = 1+mod(iJN,nstim);
%             cross_vectJN(iJN,fb_indx) = (sum(durs)*fstim_spike(fb_indx,iff) - sum(durs([1:(iJN-1) (iJN+1):end]))*stim_spike_JNf(fb_indx,iff,iJN))/durs(iJN)...
%                 - sum(durs(jn2_index))*stim_spike_JNf(fb_indx,iff,n_jn2))/durs(n_jn2);

        end
    end

    % do an svd decomposition
    ranktest(1,iff)=rank(stim_mat,ranktol);
    %[u,s,v] = svd(stim_mat);
    [u,s,v] = do_cached_calc('svd',stim_mat);
    tots = s(1,1);
    cums(iff,2) = s(1,1);
    for ii=2:nb
        tots = tots + s(ii,ii);
        cums(iff,ii+1) = cums(iff,ii) + s(ii,ii);
    end
    is = zeros(nb,nb);
    for ii=1:nb+1
        cums(iff,ii) = cums(iff,ii)/tots;
    end

    % ncutt_off = round(1.0/cums(iff,2))+1;
    % This is an adhoc way
    for ii=1:nb
        if ii>ranktest(1,iff)
            is(ii,ii)=(1.0/ranktol)*exp(-(ii-ranktest(1,iff))^2/8);
        else
            is(ii,ii)=1.0/s(ii,ii);
        end
    end
    
    % This is ridge regression
     for ii=1:nb
             is(ii,ii) = 1.0/(s(ii,ii) + ranktol);
     end
     fprintf(1,'Ridge regression\n');

    h = v*is*(u'*cross_vect);
    for  iJN=1:nJN
        hJN(iJN,:) = (v*is*(u'*cross_vectJN(iJN,:).')).';
    end

    for ii=1:nb
        ffor(ii,iff) = h(ii);
        fforJN(ii,iff,:) = hJN(:,ii);
        if iff ~= 1
            ffor(ii,nt+2-iff) = conj(h(ii));
            fforJN(ii,nt+2-iff,:) = conj(hJN(:,ii));
        end
    end
end

% FET: 8/10/2005 Make a window to enforce causality
nt2 = (nt-1)/2;
xval = -nt2:nt2;
wcausal = (atan(xval)+pi/2)/pi;

for ii=1:nb
    strfH(ii,:) = real(ifft(ffor(ii,:))).*wcausal;
    for iJN=1:nJN
        strfHJN(ii,:,iJN) = real(ifft(fforJN(ii,:,iJN))).*wcausal;
    end
end

% strfHJN_mean = mean(strfHJN,3);
% 7/12/05
strfHJN_var = zeros(size(strfHJN));
%strfHJN_nJN = zeros(size(strfHJN_mean));

% =======================================================
%  JXZ: 7/12/2005
%  Recalculate strfHJN_std
% =======================================================
% strfHJN_std = zeros(size(strfHJN));

% JXZ, 10/4/2005
% Remove divided by zero warning for nJN = 1
% Mar 15th 2011 FET back to the definition of JN std...

if nJN > 1
    for iJN=1:nJN
        strfHJN_var = strfHJN_var + (squeeze(strfHJN(:,:,iJN)) - strfH).^2;
        % strfHJN_std(:,:,iJN) = std(strfHJN(:,:,[1:(iJN-1) (iJN+1):nJN]),0,3)/sqrt(nJN-1);
    end
end
strfHJN_var = (1-1/nJN)*strfHJN_var;
strfHJN_std = squeeze(strfHJN_var.^.5);

% =======================================================
% Save the result into Output/files
% =======================================================
if save_flag == 1
    currentPath = pwd;

    outputPath = DF_PARAMS.outputPath;
    if ~isempty(outputPath)
        cd (outputPath);
    else
        disp('Saving output to Output Dir.');
        stat = mkdir('Output');
        cd('Output');
        outputPath = pwd;
    end

    save('strfH.mat', 'strfH');
    save('strfH_std.mat', 'strfHJN_std');
    for iJN=1:nJN
        filename = sprintf('strfHJN%d.mat',iJN);
        strfHJN_nJN = strfHJN(:,:,iJN);
        save(filename, 'strfHJN_nJN');
    end
    cd(currentPath);
end

% =======================================================
% END OF CAL_STRF
% =======================================================

