function [strfH, strfHJN, strfHJN_std] = df_cal_StrfSep(CSspace, CStime, stim_spike,...
    stim_spike_JN, tol, save_flag)
%
%
%  created SVD 7/28/03 ripped off of cal_Strf.m
%
% function [strfH, strfHJN, strfHJN_std] = cal_StrfSep(CSspace, CStime, stim_spike,...
%                                               stim_spike_JN, tol, save_flag)
%
%      -- Calculate strf, JackKnifed strf and contour JN strf
%     Input:
%         CSspace: stim spatial auto-correlation
%         CStime: stim temporal auto-correlation
%         stim_spike: stim_spike cross-correlation
%         stim_spike_JN: JackKnifed stim_spike cross-correlation
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

global DF_PARAMS

% =======================================================
% Check if we have all input
% =======================================================
if ~exist('save_flag','var')
    save_flag = 0;
end


% ========================================================
% Find the maximum norm of all the matrices
% ========================================================

% get sizes of various dimensions (time,space,jackknife)
% stim_spike is nt X nb (nbands)
% stim_spike_JN is nt X nb (nbands) X nJN
stim_size=size(CSspace);
stim_spike_size=size(stim_spike);
stim_spike_JNsize=size(stim_spike_JN);

nb=stim_spike_JNsize(1);
nt=stim_spike_JNsize(2);
nJN=size(stim_spike_JN,3);

% initialize to raw stim-resp cross correlation
strfH=stim_spike;
strfHJN=stim_spike_JN;

% do temporal normalization

% create matrix out of temporal autocorrelation
% assumes that CStime is a 1 X 2*nt-1 vector
if length(CStime)<2*nt-1,
    disp('error: CStime must be >= 2*nt-1 in length');
    return
end

tCStime=zeros(nt);
for uu=1:nt,
    tCStime(uu,:)=CStime((nt+1-uu):(nt*2-uu));
end

stimnorm=norm(tCStime);
ranktol=tol*stimnorm;
ranktest=rank(tCStime,ranktol);
alien_space_file = DF_PARAMS.alien_space_file;
use_alien_space = DF_PARAMS.use_alien_space;
global the_checksum

outputPath = DF_PARAMS.outputPath;
Tol_val = DF_PARAMS.Tol_val;
if use_alien_space
    cached_dir = dir_of_caches;
    loaded_alien_pointer = load(alien_space_file);
    if ~isfield(loaded_alien_pointer,'original_time_subspace_checksum')
        msg = ['Error: STRFPAK expected the file "' alien_space_file '"' char(10) 'to be a well-formed subspace file.'];
        errordlg(msg);
        error(msg);
    end
    hash_of_usv = loaded_alien_pointer.original_time_subspace_checksum;
    to_load_file = fullfile(cached_dir,[hash_of_usv,'.mat']);
    if exist(to_load_file)
        loaded = load(to_load_file);
        u = loaded.out1;
        s = loaded.out2;
        v = loaded.out3;
        clear loaded;
    else
        msg = ['Error: STRFPAK needed a cache file of the alien subspace to use, ' char(10) ...
            'but it has been deleted since the beginning of this STRFPAK session.' char(10) ...
            'Try increasing the size of your cache to fix this problem.'];
        errordlg(msg);
        error(msg);
    end
else
    [u,s,v] = df_do_cached_calc('svd',tCStime);
    original_time_subspace_checksum = the_checksum;
    original_subspace_tol_vals = Tol_val;
    save(fullfile(outputPath,'subspace.mat'),'original_subspace_tol_vals','original_time_subspace_checksum');
end

tots = s(1,1);
cums(2) = s(1,1);
for i=2:nt
    tots = tots + s(i,i);
    cums(i+1) = cums(i) + s(i,i);
end
cums = cums./tots;
is = zeros(nt,nt);

% ncutt_off = round(1.0/cums(2))+1;
for i=1:nt
    if i>ranktest
        is(i,i)=(1.0/ranktol)*exp(-(i-ranktest)^2/8);
    else
        is(i,i)=1.0/s(i,i);
    end
end

strfH = (v*is*(u'*strfH.')).';
for  iJN=1:nJN
    strfHJN(:,:,iJN) = (v*is*(u'*strfHJN(:,:,iJN).')).';
end

% do spatial normalization... analogous to spatial, but across the
% other dimension. the only real difference is replacing tCStime
% with CSspace and not transposing strfH before the decorrelation

stimnorm=norm(CSspace);
ranktol=tol*stimnorm;
ranktest=rank(CSspace,ranktol);
if use_alien_space
    cached_dir = dir_of_caches;
    loaded_alien_pointer = load(alien_space_file);
    if ~isfield(loaded_alien_pointer,'original_subspace_checksum')
        msg = ['Error: STRFPAK expected the file "' alien_space_file '"' char(10) 'to be a well-formed subspace file.'];
        errordlg(msg);
        error(msg);
    end
    hash_of_usv = loaded_alien_pointer.original_subspace_checksum;
    to_load_file = fullfile(cached_dir,[hash_of_usv,'.mat']);
    if exist(to_load_file)
        loaded = load(to_load_file);
        u = loaded.out1;
        s = loaded.out2;
        v = loaded.out3;
        clear loaded;
    else
        msg = ['Error: STRFPAK needed a cache file of the alien subspace to use, ' char(10) ...
            'but it has been deleted since the beginning of this STRFPAK session.' char(10) ...
            'Try increasing the size of your cache to fix this problem.'];
        errordlg(msg);
        error(msg);
    end
else
    [u,s,v] = df_do_cached_calc('svd',CSspace);
    original_subspace_checksum = the_checksum;
    save(fullfile(outputPath,'subspace.mat'),'original_subspace_checksum','-APPEND');

end

tots = s(1,1);
cums(2) = s(1,1);
for i=2:nb
    tots = tots + s(i,i);
    cums(i+1) = cums(i) + s(i,i);
end
is = zeros(nb,nb);
for i=1:nb+1
    cums(i) = cums(i)/tots;
end

% ncutt_off = round(1.0/cums(2))+1;
for i=1:nb
    if i>ranktest
        is(i,i)=(1.0/ranktol)*exp(-(i-ranktest)^2/8);
    else
        is(i,i)=1.0/s(i,i);
    end
end
% FET: 8/10/2005 Make a window to enforce causality
nt2 = (nt-1)/2;
xval = -nt2:nt2;
wcausal = (atan(xval)+pi/2)/pi;

strfH = (v*is*(u'*strfH));
for  iJN=1:nJN
    strfHJN(:,:,iJN) = (v*is*(u'*strfHJN(:,:,iJN)));
    for ii = 1:nb
        strfHJN(ii,:,iJN) = strfHJN(ii,:,iJN) .* wcausal;
    end
end


% do jackknife stats (mean and stderr)
strfHJN_mean = mean(strfHJN,3);
strfHJN_var = zeros(size(strfHJN_mean));
% strfHJN_nJN = zeros(size(strfHJN_mean));


% Junli: 9/9/2005
%  recalculating strfHJN_std for 2-D data
% strfHJN_std = zeros(size(strfHJN));
% The commented code was form Junli and does not look correct to me... FET
% 2011

if nJN > 1
    for iJN=1:nJN
        nstim = nJN;
        % jn2_index = 1+mod((1:(nstim-1))+iJN,nstim);
        strfHJN_var = strfHJN_var + (strfHJN(:,:,iJN) - strfHJN_mean).^2;
        % strfHJN_std(:,:,iJN) = std(strfHJN(:,:,jn2_index),0,3)/sqrt(nJN-1);

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
