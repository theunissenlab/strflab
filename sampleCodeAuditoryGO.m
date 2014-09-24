%Typical vision setup (Not part of strflab)
%---------------------------------------

addpath(genpath('/auto/k2/share/strflabGOLD'));  % add new STRF lab path
load /auto/k2/share/strflabGO/fakedata/sampleDataAuditory.mat
% 
% stimEst=x(1:900,:);
% respEst=y(1:900);
% stimVal=x(901:1000,:);
% respVal=y(901:1000);

 x = x - mean(x,2)*ones(1,size(x,2));
 x = x ./ (std(x')'*ones(1,size(x,2)));
% stimEst = (x(:,find(assign<19))');
% respEst = mean(t(:,find(assign<19)))';
% stimVal = (x(:,find(assign > 18))');
% respVal= mean(t(:,find(assign > 18)))';
% Declaring global variables.
%-----------------------------------
global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(x',t',assign);

options = trnDirectFit;
options.time_lag_max = 100;
options.cacheDir = '/auto/k2/share/strflabGOLD/tmp';
strf = 1;
datIdx =find(assign ~=12);
[strfTrained,options]=strfOpt(strf,datIdx,options);
