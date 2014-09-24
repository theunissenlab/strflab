function [spike_time, num_trials] = preprocSpikeTime(filename,stim_dur)
%
% function [my_spiketrian] = read_spikeTime(filename,stim_dur)
%  
%    1. Convert spike arriving time to on/off flag
%    2. Subgroup them into two groups

% Created by Junli, 1/22/2004
%
% Updated by Junli, 9/22/2004
%
%    Bug 1: when subgrouping multi-trial neuron spike, forgot averaging 
%           subgroup. Fixed on 9/22/04.
%
% Modified by Junli, 6/28/2005
%  
%    Add user's selection
% Modified by Junli, 8/25/2005
%    In order to save space, we will save spike arrival time to file.
%

% Read spike arrivial time from file
fid = fopen(filename, 'r');

% % Group spike to two subgroup
% spike_time = zeros(2,stim_dur);
% 
% % Assign glaobal ntrial_proper
global num_trials
trialNum = 1;
% evenTrial = 0;
% oddTrial = 0;

line=fgetl(fid);
while ( ischar(line) )
    % Read one line from a file
    templine = sscanf(line, '%g');
    
    % Need eliminate nan invalid time points
    goodindx = find(~isnan(templine));
    templine = templine(goodindx);
    
    % Also need eliminate negative time point
    goodindx = find(templine>=0);
    templine = templine(goodindx);
    
    % Round it to right spike position
    spikeIndex = (round(templine) +1)';
    
    spike_time(trialNum,:) = zeros([1 stim_dur]);
    spikeIndex = spikeIndex(find(spikeIndex < stim_dur));
    spike_time(trialNum,spikeIndex) = 1;
    
%     % Check if spikeIndex is in stim_dur range
%     spikeIndex = spikeIndex(find(spikeIndex < stim_dur));
%     
%     %spiketrian(trialNum,:) = zeros(1,stim_dur);
%     %spiketrian(trialNum, spikeIndex) = ones(1,length(spikeIndex));
%     
%     spiketrian = zeros(1,stim_dur);
%     spiketrian(spikeIndex) = ones(1,length(spikeIndex));
%     
%     if mod(trialNum, 2) == 1
%          evenTrial = evenTrial +1;
%          my_spiketrian(1,:) = my_spiketrian(1, :) + spiketrian;
%      else
%          oddTrial = oddTrial +1;
%          my_spiketrian(2,:) = my_spiketrian(2, :) + spiketrian;
%      end
%     
    trialNum = trialNum +1;
    line=fgetl(fid);
    
end

% 9/22/04: Averaging subgroups
% my_spiketrian(1,:) = my_spiketrian(1,:)/evenTrial;
% my_spiketrian(2,:) = my_spiketrian(2,:)/oddTrial;
num_trials = trialNum -1;

