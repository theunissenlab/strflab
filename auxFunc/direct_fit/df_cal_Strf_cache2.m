function [strfH, strfHJN, strfHJN_std] = df_cal_Strf_cache2(fstim, fstim_spike,...
    stim_spike_JNf, stim_size, stim_spike_size, stim_spike_JNsize,...
    nb, nt, nJN, tol, big_u,big_s,big_v,max_s,big_u_alien,big_s_alien,big_v_alien,max_s_alien,save_flag)
%
%  [strfH, strfHJN, strfHJN_std] = df_cal_Strf_cache2(fstim, fstim_spike,...
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
% Modified by FET
% 9/30/2010
% 1. Made the regularization the one given by ridge regression
% 1b. still need to do this for the "alien space" routine
% 2. I changed back the JN becuase it is already in delete one space...

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
% NOT USED: alien_min_tol = 1e-6; % The tol value to use on the current space in conjunction with the tol values in the alien space.  Needed for numerical reasons.

% ========================================================
% Find the maximum norm of all the matrices
% ========================================================
for iff=1:nf
    stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
    stim_mat(find(tril(test)))=conj(fstim(:,iff));
    stim_mat=stim_mat-diag(diag(stim_mat))+(stim_mat');

    %     nc = 1;
    %     for fb_indx=1:nb
    %         for fb_indx2=fb_indx:nb
    %             stim_mat(fb_indx,fb_indx2) = fstim(nc,iff);
    %             if fb_indx ~= fb_indx2
    %                 stim_mat(fb_indx2,fb_indx) = conj(fstim(nc,iff));
    %             end
    %             nc = nc +1;
    %         end
    %     end
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
            cross_vectJN(iJN,fb_indx) = stim_spike_JNf(fb_indx,iff,iJN);  % FET 09/2010 - the cross-correlation is already normalized?
            %Normalized by the duration of each individual stim
            % cross_vectJN(iJN,fb_indx) = (sum(durs)*fstim_spike(fb_indx,iff) - sum(durs([1:(iJN-1) (iJN+1):end]))*stim_spike_JNf(fb_indx,iff,iJN))/durs(iJN);

        end
    end
    use_alien_space = DF_PARAMS.use_alien_space;
    if use_alien_space
        % first project into the alien space; do not do an inverse.
        u_a = big_u_alien{iff};
        s_a = big_s_alien{iff};
        v_a = big_v_alien{iff};
        is_a = zeros(nb,nb);
        ranktest_a(1,iff)=sum(diag(s_a)>tol*max_s_alien);
        % ncutt_off = round(1.0/cums(iff,2))+1;
        for ii=1:nb  % Just project; don't normalize/divide.
            if ii>ranktest_a(1,iff)
                is_a(ii,ii)=(1.0)*exp(-(ii-ranktest_a(1,iff))^2/8);
            else
                is_a(ii,ii)=1.0;
            end
        end
        %cross_vect = v*is*(u'*cross_vect);
        % Now do the matrix inversion
        ranktest(1,iff)=rank(stim_mat);  %Still useful, since we don't want division by exactly 0.  We'll let MATLAB handle numerical issues and machine precision.
        u = big_u{iff};
        s = big_s{iff};
        v = big_v{iff};
        is = zeros(nb,nb);
        for ii=1:nb  % Just divide; don't project; this is accomplished by making ranktest very large as above.
            if ii>ranktest(1,iff)
                is(ii,ii)=(1.0/s(ranktest(1,iff),ranktest(1,iff)))*exp(-(ii-ranktest(1,iff))^2/8);
            else
                is(ii,ii)=1.0/s(ii,ii);
            end
        end

        h = v_a*is_a*(u_a'*(v*is*(u'*cross_vect)));
        for  iJN=1:nJN
            hJN(iJN,:) = (v_a*is_a*(u_a'*(v*is*((u'*cross_vectJN(iJN,:).'))))).';
        end

        for ii=1:nb
            ffor(ii,iff) = h(ii);
            fforJN(ii,iff,:) = hJN(:,ii);
            if iff ~= 1
                ffor(ii,nt+2-iff) = conj(h(ii));
                fforJN(ii,nt+2-iff,:) = conj(hJN(:,ii));
            end
        end

    else
        if 0 % debugging Feb 5 2007
            % first project into the alien space; do not do an inverse.
            u_a = big_u{iff};
            s_a = big_s{iff};
            v_a = big_v{iff};
            is_a = zeros(nb,nb);
            ranktest_a(1,iff)=rank(stim_mat,ranktol);

            % ncutt_off = round(1.0/cums(iff,2))+1;
            for ii=1:nb  % Just project; don't normalize/divide.
                if ii>ranktest_a(1,iff)
                    is_a(ii,ii)=(1.0)*exp(-(ii-ranktest(1,iff))^2/8);
                else
                    is_a(ii,ii)=1.0;
                end
            end

            %cross_vect = v*is*(u'*cross_vect);
            % Now do the matrix inversion
            ranktest(1,iff)=rank(stim_mat,max_s*alien_min_tol);  %Still useful, since we don't want division by exactly 0.  We'll let MATLAB handle numerical issues and machine precision.
            u = big_u{iff};
            s = big_s{iff};
            v = big_v{iff};
            is = zeros(nb,nb);

            % ncutt_off = round(1.0/cums(iff,2))+1;
            for ii=1:nb  % Just divide; don't project; this is accomplished by making ranktest very large as above.
                if ii>ranktest(1,iff)
                    is(ii,ii)=(1.0/s(ranktest(1,iff),ranktest(1,iff)))*exp(-(ii-ranktest(1,iff))^2/8);
                else
                    is(ii,ii)=1.0/s(ii,ii);
                end
            end

            h = v_a*s_a*(u_a'*(v*is*(u'*cross_vect)));
            for  iJN=1:nJN
                try
                hJN(iJN,:) = (v_a*s_a*(u_a'*(v*is*((u'*cross_vectJN(iJN,:).'))))).';
                catch
                    disp(['size of 1st bracket:' size((u'*cross_vectJN(iJN,:).'))]);
                    error(lasterr);
                end
            end

            for ii=1:nb
                ffor(ii,iff) = h(ii);
                fforJN(ii,iff,:) = hJN(:,ii);
                if iff ~= 1
                    ffor(ii,nt+2-iff) = conj(h(ii));
                    fforJN(ii,nt+2-iff,:) = conj(hJN(:,ii));
                end
            end

        else
            % do an svd decomposition
            ranktest(1,iff)=rank(stim_mat,ranktol);
            %[u,s,v] = svd(stim_mat);
            %[u,s,v] = do_cached_calc('svd',stim_mat);
            u = big_u{iff};
            s = big_s{iff};
            v = big_v{iff};
            is = zeros(nb,nb);

            % ncutt_off = round(1.0/cums(iff,2))+1;
%             for ii=1:nb
%                 if ii>ranktest(1,iff)
%                     is(ii,ii)=(1.0/ranktol)*exp(-(ii-ranktest(1,iff))^2/8);
%                 else
%                     is(ii,ii)=1.0/s(ii,ii);
%                 end
%             end
            % This is ridge regression
            for ii=1:nb
                is(ii,ii) = 1.0/(s(ii,ii) + ranktol);
            end
            % fprintf(1,'Ridge regression\n');

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
    end
end

% FET: 8/10/2005 Make a window to enforce causality
nt2 = (nt-1)/2;
xval = -nt2:nt2;
global strictcausal
if strcmp(strictcausal,'true')
    wcausal = xval>0;
else
    wcausal = (atan(xval)+pi/2)/pi;
end

for ii=1:nb
    strfH(ii,:) = real(ifft(ffor(ii,:))).*wcausal;
    for iJN=1:nJN
        strfHJN(ii,:,iJN) = real(ifft(fforJN(ii,:,iJN))).*wcausal;
    end
end

% strfHJN_mean = mean(strfHJN,3);
% 7/12/05
strfHJN_var = zeros(size(strfH));
%strfHJN_nJN = zeros(size(strfHJN_mean));

% =======================================================
%  JXZ: 7/12/2005
%  Recalculate strfHJN_std
% =======================================================
% strfHJN_std = zeros(size(strfHJN));

% JXZ, 10/4/2005
% Remove divided by zero warning for nJN = 1
if nJN > 1
    for iJN=1:nJN
        strfHJN_var = strfHJN_var + (squeeze(strfHJN(:,:,iJN)) - strfH).^2;
       % strfHJN_std(:,:,iJN) = std(strfHJN(:,:,[1:(iJN-1) (iJN+1):nJN]),0,3);  % This is already in standard errors
    end
end
strfHJN_var = (1-1/nJN)*strfHJN_var;  % This is the correct definition for the standard error of the mean
strfHJN_std = strfHJN_var.^.5;

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

