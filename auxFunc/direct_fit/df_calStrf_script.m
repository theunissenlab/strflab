% =========================================
% matlab script file for calStrf
% ===========================================
% Load up all files and Display for checking
% still need to look at .rec files for checking consistency
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
%

global DF_PARAMS

nstd_val = 0.5;

% ===========================================
% FFT Auto-correlation and Cross-correlation
% ===========================================

%pack;
checksum_fft_ACC= df_checksum(checksum_CrossCorr,round(TimeLag*ampsamprate/1000),nstd_val);
TimeLagUnit = DF_PARAMS.TimeLagUnit;

if strcmp(TimeLagUnit, 'msec')
    twindow = round(TimeLag*ampsamprate/1000);
elseif strcmp(TimeLagUnit, 'frame')
    twindow = round(TimeLag);
end
nt = 2*twindow + 1;

if cache_crosscorr
    [fstim, fstim_spike, stim_spike_JNf] = df_do_locally_cached_calc_checksum_known(...
        df_get_local_cache_dir,'df_fft_AutoCrossCorr',checksum_fft_ACC,CS,...
        CSR,CSR_JN, twindow, NBAND, nstd_val);
else
    [fstim, fstim_spike, stim_spike_JNf] = df_fft_AutoCrossCorr(CS,...
        CSR,CSR_JN, twindow, NBAND, nstd_val);
end
disp('Done df_fft_AutoCrossCorr.');

% clear some memory
clear CSR CSR_JN CS

%pack;
% ===========================================
%  Prepare for call STRF_calculation
% ===========================================
nb = NBAND;
nJN = length(DS);
stim_size = size(fstim);
stim_spike_size = size(fstim_spike);
stim_spike_JNsize = size(stim_spike_JNf);


% ===========================================
% Get tolerance values
% ===========================================
Tol_val = DF_PARAMS.Tol_val;
ntols = length(Tol_val);
outputPath = DF_PARAMS.outputPath;
if isempty(outputPath)

    disp('Saving output to Output Dir.');
    stat = mkdir('Output');
    outputPath = fullfile(pwd, 'Output');
end

% ===========================================
fprintf('Calculating STRF for each tol value...\n');
nf = (nt-1)/2 +1;
hashes_of_stims = df_create_stim_cache_file(outputPath,DS);
hashes_of_spikes = df_create_spike_cache_file(outputPath,DS);

use_more_memory = 1;  %Turning this off will break the re-using subspace special option.
use_alien_space = DF_PARAMS.use_alien_space;
alien_space_file = DF_PARAMS.alien_space_file;
if length(use_alien_space) == 0
    use_alien_space = 0;
end
if use_alien_space
    cached_dir = df_dir_of_caches;
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
        big_u_alien = loaded.out1;
        big_s_alien = loaded.out2;
        big_v_alien = loaded.out3;
        max_stimnorm_alien = loaded.out4;
        big_stim_mat_alien = loaded.out5;
        max_s_alien = loaded.out6;
        clear loaded;

        if length(big_u_alien) ~= nf
            msg = ['Error: mismatch in time domain between the alien subspace and the current subspace.'];
            errordlg(msg)
            error(msg);
        end
        if size(big_u_alien{1},1) ~= nb
            msg = ['Error: mismatch in space domain between the alien subspace and the current subspace.'];
            errordlg(msg)
            error(msg);
        end
    else
        msg = ['Error: STRFPAK needed a cache file of the alien subspace to use, ' char(10) ...
            'but it has been deleted since the beginning of this STRFPAK session.' char(10) ...
            'Try increasing the size of your cache to fix this problem.'];
        errordlg(msg);
        error(msg);
    end
end
if use_more_memory
    big_u = {};
    big_s = {};
    big_v = {};

    %checksum_usv = df_checksum(df_load_function_text('df_make_big_usv'),hashes_of_spikes,hashes_of_stims,twindow,round(TimeLag*ampsamprate/1000),nstd_val);
    checksum_usv = df_checksum(df_load_function_text('df_make_big_usv'),hashes_of_stims,twindow,round(TimeLag*ampsamprate/1000));
    if 1  %This is being forced so that re-using subspaces will not work only in gui mode.
        [big_u,big_s,big_v,max_stimnorm,big_stim_mat,max_s] = df_do_cached_calc_checksum_known('df_make_big_usv',checksum_usv,nb,nf,fstim); %max_s needed here or it won't be cached
        if ~use_alien_space
            original_subspace_tol_vals = Tol_val;
            original_subspace_checksum = checksum_usv;
            save(fullfile(outputPath,'subspace.mat'),'original_subspace_tol_vals','original_subspace_checksum');
        end
    else
        [big_u,big_s,big_v,max_stimnorm,big_stim_mat,max_s] = df_make_big_usv(nb,nf,fstim); 
    end
end

%if cache_crosscorr
checksum_cal_strf = df_checksum(checksum_fft_ACC,stim_size,stim_spike_size,stim_spike_JNsize, nb, nt, nJN);
if use_alien_space
    checksum_cal_strf = df_checksum(checksum_cal_strf,loaded_alien_pointer);
end
%end
for itol=1:ntols
    tol=Tol_val(itol);
    if use_alien_space
        if ~any(loaded_alien_pointer.original_subspace_tol_vals == tol)
            msg = ['Error: STRFPAK was asked to calculate a STRF using the tol value ' num2str(tol) char(10) ...
                'using a previously-computed subspace, but that subspace used only the tol values' char(10) ...
                num2str(loaded_alien_pointer.original_subspace_tol_vals) '.' char(10) char(10) ...
                '(If you intended to use a new tol value with the alien subspace, disable this alarm by typing' char(10) ...
                '"edit ' mfilename '" and commenting out this error.  STRFPAK won''t balk, but' char(10) ...
                'we can''t think of a reason to use an alien subspace other than to get the same regularization bias' char(10)...
                'for two different stim ensembles, and using different tol values defeats the point here.)'];
            errordlg(msg);
            error(msg);
        end
    end

    fprintf('Now calculating STRF for tol_value: %g\n', tol);

    % =======================================
    % Calculate strf for each tol val.
    % =======================================
    checksum_this_tol = df_checksum(checksum_cal_strf,df_load_function_text('df_cal_Strf'),df_load_function_text('df_cal_Strf_cache2'),tol);
    if (~ cache_crosscorr)
        if ~use_more_memory
            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_cal_Strf(fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol);
        else
            %         [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_cal_Strf_use_cache(fstim,...
            %             fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
            %             stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_stimnorm,big_stim_mat);
            if ~ use_alien_space
            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_cal_Strf_cache2(fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_s);
            else
                            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_cal_Strf_cache2(fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_s,big_u_alien,big_s_alien,big_v_alien,max_s_alien);
            end

        end
    else

        if ~use_more_memory
            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_Strf',checksum_this_tol,fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol);
        else
            %         [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_cal_Strf_use_cache(fstim,...
            %             fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
            %             stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_stimnorm,big_stim_mat);
            if ~ use_alien_space
            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_Strf_cache2',checksum_this_tol,fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_s);
            else
                            [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_Strf_cache2',checksum_this_tol,fstim,...
                fstim_spike, stim_spike_JNf,stim_size, stim_spike_size,...
                stim_spike_JNsize, nb, nt, nJN, tol,big_u,big_s,big_v,max_s,big_u_alien,big_s_alien,big_v_alien,max_s_alien);

            end
        end
    end


    fprintf('Done calculation of STRF for tol_value: %g\n', tol);


    sfilename = sprintf('strfResult_Tol%d.mat',itol);
    strfFiles = fullfile(outputPath,sfilename);
    save(strfFiles, 'STRF_Cell', 'STRFJN_Cell', 'STRFJNstd_Cell');
    if cache_crosscorr  %  You could get a checksum directly from the STRF, but this is much faster and just as unique.
        strf_checksum = checksum_this_tol;
        posslash = findstr(strfFiles,filesep);
        the_dir = strfFiles(1:(posslash(end)));
        the_name = strfFiles((posslash(end)+1):end);
        strf_hash_filename = [the_dir 'hash_of_' the_name];
        save(strf_hash_filename,'strf_checksum');
    end
    clear STRF_Cell STRFJN_Cell STRFJNstd_Cell
end

% ===========================================
%  END of df_calStrf_script.m
% ===========================================

