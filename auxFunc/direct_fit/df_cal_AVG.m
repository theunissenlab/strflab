function [stim_avg, Avg_psth_out,psth,constmeanrate,errFlg] = df_cal_AVG(DDS, nband, psth_option, lin_flag, sil_window)
%
% [stim_avg, Avg_psth, psth] = df_cal_AVG(DDS, lin_flag, sil_window)
%     -- Calculate the average stimulus value across time.
%     -- Calculate average of psth over all trials 
%     -- Calculate the psth of response file over all trials.
%   Input:
%      DDS(required input):  the cell of each data struct 
%                               that contains four fields:
%               stimfiles  - stimulus file name
%               respfiles  - response file name
%               nlen    - length of time domain
%               ntrials    - num of trials
%              e.g. DDS{1} = struct('stimfiles', 'stim1.dat', 'respfiles', 
%                    'resp1.dat', 'nlen', 1723, 'ntrials', 20);
%      lin_flag: the flag to show whether we need take log on data
%              e.g. lin_flag = 0 if we need take log, 1(default) if otherwise 
%      sil_window: the interval for not counting when preprocessing data
%              e.g. sil_window = 0(default) 
%   Output:
%      stim_avg:  average stimulus value which is only function of space.
%      Avg_psth: average psth 
%      psth: the cell of psth for one data pair which is only func of time 
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

% created by JXZ
% Sep. 25, 2002
%
%
% Junli: 6/30/2005
%    Add parameter 'psth_option' to specifiy psth noise removal option

global DF_PARAMS

% ========================================================
% check whether we have valid required input
% ========================================================
errFlg = 0;
if isempty(DDS) 
    fprintf('ERROR: Please enter non-empty data filename');
    errFlg = 1;
    return
end

NBAND = DF_PARAMS.NBAND;
if ~exist('nband')
    if isempty(NBAND)
        fprintf('You need assign variable NBAND first.');
        errFlg = 1;    
        return
    end
    nband = NBAND;
end


% ========================================================
% Junli: 6/30/2005
%    Add parameter 'psth_option' to specifiy psth noise removal option
% ========================================================
timevary_PSTH = DF_PARAMS.timevary_PSTH;
if ~exist('psth_option')
    if timevary_PSTH == 0
         psth_option = 0;
    else 
        psth_option = 1;
    end
end

% ========================================================
% if user does not provide other input, we set default.
% ========================================================
% lin_flag refers to linear data if lin_flag = 1.
% otherwise, we take logrithm on data if lin_flag = 0.
if ~exist('lin_flag')
    lin_flag = 1;
end

if ~exist('sil_window')
    sil_window = 0;
end

% ========================================================
% initialize output and declare local variables 
% ========================================================
% find out the total number of data files
ampsamprate = DF_PARAMS.ampsamprate;
respsamprate = DF_PARAMS.respsamprate;
ndata_files = length(DDS);

stim_avg = zeros(nband, 1);
count_avg = 0;
tot_trials = 0;
psth = {};
timevary_psth = [];
Avg_psth = 0;

% ========================================================
% calculate the output over all the data file 
% ========================================================
for n = 1:ndata_files
    % load stimulus files
    stim_env = df_Check_And_Load(DDS{n}.stimfiles);
    this_len = size(stim_env,2);  % get stim duration
    
    % load response files
    rawResp = df_Check_And_Load(DDS{n}.respfiles);
   
    if iscell(rawResp)
        spiketrain = zeros(DDS{n}.ntrials,this_len);
        for trial_ind =1:DDS{n}.ntrials
            spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
        end
        if ~isempty(ampsamprate) & ~isempty(respsamprate)
            newpsth = resample(spiketrain', ampsamprate, respsamprate);
        else
            newpsth = spiketrain;
        end
        newpsth = newpsth'; % make sure new response data is trials x T.
        newpsth(find(newpsth < 0)) = 0;
        psth_rec = newpsth;
    else
        psth_rec = rawResp;
    end
    nt = min(size(stim_env, 2), size(psth_rec,2));
    % take logrithm of data based on lin_flag 
    if lin_flag == 0
        stim_env = log10(stim_env + 1.0);
    end
    
    % calculate stim avg
    %
    % Before Do calculation, we want to check if we got the correct input
    tempXsize = size(stim_env,1);
    if tempXsize ~= nband
        fprintf('Data error, first data file needs to be stimuli, second needs to be response.');
        errFlg = 1;
        return
    end
    
    stim_avg = stim_avg + sum(stim_env*DDS{n}.ntrials, 2);
    count_avg = count_avg +(nt + 2*sil_window)*DDS{n}.ntrials;  
    
    % then calculate response_avg
    if DDS{n}.ntrials > 1        
        temp= mean(psth_rec);
        psth{n} = temp(1:nt);
    else
        
        psth{n} = psth_rec(1:nt);
    end
    
    tot_trials = tot_trials + nt + sil_window;
    
    % calculate the total spike/response avg.
    Avg_psth = Avg_psth + sum(psth{n}(1:nt));
    
    timevary_psth = [timevary_psth size(psth{n},2)];
    
    % clear workspace
    clear stim_env
    clear psth_rec
end

% Junli: 6/29/2005
% -------------------------------------------------------
%  Calculating Time Varying mean firing rate 
% -------------------------------------------------------
max_psth_indx = max(timevary_psth);
whole_psth = zeros(length(psth), max_psth_indx);
count_psth = zeros(length(psth), max_psth_indx);

for nn = 1:length(psth) 
    whole_psth(nn,:) = [psth{nn}.*DDS{nn}.ntrials zeros(1,max_psth_indx-size(psth{nn},2))];
    count_psth(nn,:) = [ones(1,size(psth{nn},2)).*DDS{nn}.ntrials zeros(1,max_psth_indx-size(psth{nn},2))];      
end

sum_whole_psth = sum(whole_psth);
sum_count_psth = sum(count_psth);

% Make Delete one averages and smooth at window = 41 - This needs to be a
% parameter
smooth_rt = DF_PARAMS.smooth_rt;
if isempty(smooth_rt)
    smooth_rt = 41;
end
psthsmoothconst = smooth_rt;
if mod(psthsmoothconst, 2) == 0
    cutsize = 0;
else
    cutsize = 1;
end
halfwinsize = floor(psthsmoothconst/2);
wind1 = hanning(psthsmoothconst)/sum(hanning(psthsmoothconst));  

for nn = 1:length(psth)
    count_minus_one = sum_count_psth-count_psth(nn,:);
    count_minus_one(find(count_minus_one==0))=1;  % Set these to 1 to prevent 0/0
    whole_psth(nn,:) = (sum_whole_psth - whole_psth(nn,:))./count_minus_one;
    svagsm=conv(whole_psth(nn,:), wind1);
    whole_psth(nn,:) = svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize);
end


% ========================================================
% save the stim_avg into the data file
% ========================================================
currentPath = pwd;
outputPath = DF_PARAMS.outputPath;
if ~isempty(outputPath)
    cd (outputPath);
else
    % Take care of empty outputPath case
    currentPath = pwd;
    fprintf('No output path specified for intermediate results, defaulting to %s', currentPath);
    outputPath = currentPath;
end

stim_avg = stim_avg/count_avg;

% Smoothing time varying mean firing rate using big smooth window
ampsamprate = DF_PARAMS.ampsamprate;

constmeanrate = Avg_psth / tot_trials;
if psth_option == 0
    Avg_psth_out = constmeanrate;      
else
    Avg_psth_out = whole_psth;
end
Avg_psth = whole_psth;
save(fullfile(outputPath,'stim_avg.mat'), 'stim_avg', 'Avg_psth','constmeanrate');

% ========================================================
% Junli: 6/30/2005
%     Display for debugging
% ========================================================
%if psth_option == 1
if 0
     % Display for debugging     
     if mod(psth_smoothconst, 2) == 0
        cutsize = 0;
    else
        cutsize = 1;
    end
    %Show smoothed version of psth at window = 15
    halfwinsize = floor(psth_smoothconst/2);
    wind1 = hanning(psth_smoothconst)/sum(hanning(psth_smoothconst)); 
    svagsm=conv(mean(whole_psth,1),wind1);
    figure
    plot(svagsm(halfwinsize+cutsize:length(svagsm)-halfwinsize)*ampsamprate);
    xlabel('Time(ms)')
    ylabel('Time varying mean rate (Spikes/s)')
   
end

cd(currentPath);
% ========================================================
% END OF CAL_AVG
% ========================================================

