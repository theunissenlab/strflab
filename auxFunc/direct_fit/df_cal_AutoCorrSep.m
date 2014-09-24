function [CSspace, CStime, errFlg] = df_cal_AutoCorrSep(DS, stim_avg,twindow, nband, JN_flag)
%
%

global DF_PARAMS

% ========================================================
% check whether we have valid required input
% ========================================================
errFlg = 0;
if isempty(DS)
    fprintf('ERROR: Please enter non-empty data filename\n');
    errFlg = 1;
    return
   
end

if ~exist('JN_flag')
    JN_flag = 0;
end

% ========================================================
% check whether stim_avg has been calculated or not
% ========================================================
if isempty(stim_avg)
    % calculate stim_avg and avg_psth and psth 
    [stim_avg, avgr, psth] = df_cal_AVG(DS, nband);
end

% ========================================================
% initialize the output and set its size
% ========================================================
% get the total data files
filecount = length(DS); 

% temporal axis range
tot_corr = diff(twindow) + 1;

% spatial axis range
spa_corr = (nband * (nband - 1))/2 + nband;

% initialize CS and CSJN variable

CSspace = zeros(nband,nband);
CStime = zeros(1,tot_corr*2-1);
CSspace_ns = 0;
CStime_ns = zeros(1, tot_corr*2-1);

% ========================================================
% do calculation. The algorithm is based on FET's dcp_stim.c
% ========================================================

%Visually give calculate auto_corr status
for fidx = 1:filecount     % loop through all data files
    
    % load stimulus file
    stim_env = df_Check_And_Load(DS{fidx}.stimfiles);

    nlen = DS{fidx}.nlen;
    xb = 1;
    stimval = zeros(nband, nlen);
    
    % Check if input data are chosen properly.
    thisLength = size(stim_env, 2);
    if thisLength < nlen
        answ= questdlg(['Data Error: Please check your input data by clicking ',...
     '"Get Files" Button in the main window: The first data file need ',...
     'to be stimuli and the second data file need to be its corresponding',...
     ' response file. If you made a mistake, please type "clear all" ',...
     ' or hit "reset" button first and then choose input data again.',...
     ' Otherwise, we will truncate the longer length. Do you want to continue?'],...
     'Input Data Error','Yes','No','No');
       switch answ
       case 'No'
           errFlg = 1;
           return
       end
    end
      
    nlen = min(nlen, thisLength);
    
    % subtract mean of stim from each stim 
    for tgood = 1:nlen
        stimval(:, tgood) = stim_env(1:nband, tgood) - stim_avg(1:nband);
    end
    
    
    % do autocorrelation calculation

    % spatial autocorr - time lag=0
    % gives nband X nband matrix
    CSspace=CSspace + DS{fidx}.ntrials .* stimval*stimval';
    
    % Count the total trials for later normalization 
    CSspace_ns = CSspace_ns + DS{fidx}.ntrials * nlen;
    
    % done with space

    % temporal autocorr - averaged over all spatial channels
    
    for ib1 = 1:nband

       % NEW version of algorithm by using xcorr 
       CStime = CStime + xcorr(stimval(ib1,:), stimval(ib1, :),...
                               tot_corr-1) .* DS{fidx}.ntrials;
       
    end                 % END of ib1
    
    % Count the total trials for later normalization 
    lengthVec = ones(1, nlen);
    CStime_ns = CStime_ns + DS{fidx}.ntrials * nband * ...
        xcorr(lengthVec, lengthVec, tot_corr-1);
   

   % clear workspace
   clear stim_env
   clear stimval
   clear lengthVec
  
end		        % END of fidx


disp('Done auto-correlation  calculation');

% ========================================================
% To normalize CS by CS_ns: 
%   if CS_ns ~= 0
%      CS = CS /CS_ns
%   end
% ========================================================

% normalize CS matrix 
CStime=CStime./(CStime_ns + (CStime==0));
CSspace=CSspace./(CSspace_ns + (CSspace==0));

CStime=CStime./mean(diag(CSspace));

% ========================================================
% save stimulus auto-corrlation matrix into file 
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

save('Stim_autocorr.mat', 'CSspace', 'CStime'); 
cd(currentPath);
% ========================================================
% END OF CAL_AUTOCORRSEP
% ========================================================
