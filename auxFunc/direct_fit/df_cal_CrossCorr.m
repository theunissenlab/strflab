function [CSR, CSR_JN, errFlg] = df_cal_CrossCorr(DS, stim_avg, avg_psth, psth,...
                      twindow, nband, end_window, JN_flag)
%
% CSR = df_cal_CrossCorr(DS, stim_avg, avg_psth, psth, twindow,
%                            nband, end_window, JN_flag)
%     -- Calcualte stimulus spike cross-correlation.
%     Input:
%        DS: the cell of each data struct that contains four fields:
%               stimfiles  - stimulus file name
%               respfiles  - response file name
%               nlength    - length of time domain
%               ntrials    - num of trials
%              e.g. DS{1} = struct('stimfiles', 'stim1.dat', 'respfiles',
%                    'resp1.dat', 'nlen', 1723, 'ntrials', 20);
%        stim_avg: avg stimulus that used to smooth the noise
%        avg_psth: average psth over all trials and over all time
%        psth: the cell of avg. psth over trials
%        twindow: the variable to set the time interval to calculate
%                 autocorrelation. e.g. twindow=[-300 300]
%        nband: the size of spatio domain of the stimulus file
%        end_window: the time interval which dont count for data analysis 
%        JN_flag: the flag that specify whether we calculate JackKnifed CS
%                  The default value of JN_flag = 1(calulate), 0 otherwise
%   Output:
%        CSR: stimulus spike cross correlation. Its size is:
%            nband X (2*twindow +1)
%        CSR_JN: stimulus spike cross correlation. Its size is:
%            nband X (2*twindow +1)
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
% Sep. 27, 2002
% Modified by JXZ
% Jun. 27, 2005

global DF_PARAMS

% ========================================================
% check whether we have valid required input
% ========================================================
errFlg = 0;
if isempty(DS) 
    errordlg('ERROR: Please enter non-empty data filename',...
        'Input Data Error', 'modal')
    errFlg = 1;
    return
end

if ~exist('JN_flag')
    JN_flag = 1;
end

if ~exist('end_window')
    end_window = 0;
end

% ========================================================
% check whether stim_avg has been calculated or not
% ========================================================
if isempty(stim_avg) 
    % calculate avg. of stimuli and psh and total_psth
    [stim_avg, avg_psth, psth] = df_cal_AVG(DS, nband);
end

% ========================================================
% initialize the output and allocate its size
% ========================================================
% get the total data files
filecount = length(DS); 

% temporal axis range
tot_corr = diff(twindow) + 1;

% spatial axis range
spa_corr = nband;

% initialize autoCorr and autoCorrJN variable
CSR = zeros(spa_corr, tot_corr);
CSR_ns = zeros(1,tot_corr);
CSR_JN = cell(filecount, 1);
CSR_JN_ns = cell(filecount, 1);

% ========================================================
% initialize JN variables
% ========================================================
if JN_flag == 1
    for iJN=1:filecount
        CSR_JN{iJN} = zeros(spa_corr, tot_corr);
        CSR_JN_ns{iJN} = zeros(1, tot_corr);
    end
end

disp('Now doing cross-correlation calculation.');
% ========================================================
% do calculation. The algorithm is from FET's dcp_stim_spike
% ========================================================
for fidx = 1:filecount     % loop through all data files
    % load stimulus file
    stim_env = df_Check_And_Load(DS{fidx}.stimfiles);
    
    % get time length for data input set
    nlen = min(size(psth{fidx}, 2),size(stim_env, 2));
    
%     % Check if input data are chosen properly.
%     thisLength = size(stim_env, 2);
%     if thisLength < nlen
%         answ= questdlg(['Data Error: Please check your input data by clicking ',...
%                 '"Get Files" Button in the main window: The first data file need ',...
%                 'to be stimuli and the second data file need to be its corresponding',...
%                 ' response file. If you made a mistake, please type "clear all" ',...
%                 ' or hit "reset" button first and then choose input data again.',...
%                 ' Otherwise, we will truncate the longer length. Do you want to continue?'],...
%             'Input Data Error','Yes','No','No');
%         switch answ
%             case 'No'
%                 errFlg = 1;
%                 return
%         end
%     end
%     
%     nlen = min(nlen, thisLength);
    
    % subtract mean_stim from stim and mean_psth from psth
    stimval = zeros(nband, nlen);
    for tgood = 1:nlen
        stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);    
    end
    
    % 6/27/2005
    % For Time-varying firing rate
    %
    timevary_PSTH = DF_PARAMS.timevary_PSTH;
    if timevary_PSTH ==1
        
        psthval = psth{fidx}(1:nlen) - avg_psth(fidx, 1:nlen);
    else
        psthval = psth{fidx} - avg_psth;
    end
    
    % New version of algorithm for computing cross-correlation
    CSR_JN{fidx} = df_internal_cal_CrossCorr(stimval,psthval,twindow(2));
    CSR = CSR + CSR_JN{fidx};
%     for ib1 = 1:nband
%         CSR_JN{fidx}(ib1,:) = xcorr(stimval(ib1, :), psthval, twindow(2));
%         CSR(ib1, :) = CSR(ib1, :) + CSR_JN{fidx}(ib1,:);
%     end
    %CSR_JN{fidx} = CSR;
    
    % For normalization and assign the count_ns
    CSR_JN_ns{fidx} = xcorr(ones(1, nlen), ones(1, nlen), twindow(2));
    CSR_ns = CSR_ns + CSR_JN_ns{fidx};
    
    
    % clear workspace by deleting stim_env
    clear stim_env;
    clear stimval;
    clear psthval;
    
end		              % END of fidx
disp('Done calculation of cross-correlation.');

disp('Now calculating JN cross-correlation.');
% calculate JN version of cross-correlation
if JN_flag == 1

    if filecount >1
        for iJN=1:filecount

            % Count ns for each JN and normalize it later on
            CSR_JN_ns{iJN} = CSR_ns - CSR_JN_ns{iJN};
            nozero_ns = isinf( 1 ./CSR_JN_ns{iJN}) + CSR_JN_ns{iJN};

            for ib = 1:nband
                CSR_JN{iJN}(ib,:) = (CSR(ib,:) - CSR_JN{iJN}(ib,:)) ./nozero_ns;
            end

        end % End of iJN
    end % END of filecount > 1
end     % End of JN_flag

disp('Done calculation of JN cross-correlation.');


% ========================================================
% To normalize CS by CS_ns:
%   if CSR_ns ~= 0
%      CSR = CSR /CSR_ns
%   end
% ========================================================
% elminate zero in autoCorr_ns 
nozero_ns = isinf( 1 ./CSR_ns) + CSR_ns;

% normalize autoCorr matrix 
for i=1:nband
    CSR(i,:) = CSR(i,:) ./ nozero_ns;
end

% ========================================================
% save stim-spike cross correlation matrix into a file
% ========================================================
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

save('StimResp_crosscorr.mat', 'CSR');  
save('SR_crosscorrJN.mat', 'CSR_JN');
cd(currentPath);

% ========================================================
% END OF CAL_AUTOCORR
% ========================================================
 
