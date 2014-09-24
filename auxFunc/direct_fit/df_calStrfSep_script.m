% ===========================================
% matlab script file for calStrf_Sep
% ===========================================
%  It does:
%    1. Intialize variables for calStrfSep.
%    2. Load up all files
%    3. Calculate STRF and STRF_JN and STRF_JNstd for each tol value.
%    4. Save them to CELL struct that easily saved.
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
%

global DF_PARAMS

% Initialize variables for cross-correlation function cal_StrfSep
nb = NBAND;
if strcmp(TimeLagUnit, 'msec')
    nt = 2*round(TimeLag*ampsamprate/1000) +1;
else
    nt = 2*round(TimeLag) +1;
end

nJN = length(CSR_JN);

% First smooth cross-correlation by hanning window with size nt
w = hanning(nt);
stim_spike = fliplr(CSR);
for ib=1:nb
    stim_spike(ib,:)=stim_spike(ib,:).*w';
end

%  Read the JN cross-correlation
stim_spike_JN = zeros(nb, nt, nJN);
w = hanning(nt);
for iJN=1:nJN
    CSR = CSR_JN{iJN};
    for ib=1:nb
        stim_spike_JN(ib,:,iJN) = df_wrev(CSR(ib,:)) .* w';
    end
end

% Get tolerance values
Tol_val = DF_PARAMS.Tol_val;
ntols = length(Tol_val);

outputPath = DF_PARAMS.outputPath;
if isempty(outputPath)

    disp('Saving output to Output Dir.');
    stat = mkdir('Output');
    outputPath = fullfile(pwd, 'Output');
end

% =======================================
% Calculate strf for each tol val.
% =======================================
checksum_cal_StrfSep = df_checksum(checksum_CrossCorr,'df_cal_StrfSep');
use_alien_space = DF_PARAMS.use_alien_space;
alien_space_file = DF_PARAMS.alien_space_file;
if use_alien_space
    extra_to_hash = load(alien_space_file);
else
    extra_to_hash = 'Not used.';
end
for itol=1:ntols
    tol=Tol_val(itol);
    fprintf('Now calculating STRF for tol_value: %g\n', tol);

    if 0 %~cache_crosscorr
        [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = cal_StrfSep(CSspace, CStime, stim_spike,...
            stim_spike_JN, tol, 1);
    else
        the_checksum = df_checksum(checksum_cal_StrfSep,tol,df_load_function_text('df_cal_StrfSep.m'),extra_to_hash);
        [STRF_Cell, STRFJN_Cell, STRFJNstd_Cell] = df_do_locally_cached_calc_checksum_known(df_get_local_cache_dir,'df_cal_StrfSep',the_checksum,CSspace, CStime, stim_spike,...
            stim_spike_JN, tol, 1);
    end

    fprintf('Done calculation of  STRF for tol_value: %g\n', tol);

    sfilename = sprintf('strfResult_Tol%d.mat',itol);
    save(fullfile(outputPath,sfilename), 'STRF_Cell', 'STRFJN_Cell', 'STRFJNstd_Cell');

    clear STRF_Cell STRFJN_Cell STRFJNstd_Cell

end

% ====================================================================
% End of calStrfSep_script.m
% ====================================================================

