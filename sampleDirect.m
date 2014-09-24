clear global
load /auto/k5/prenger/strflab/fakedata/sampledata_nonan.mat
%rmpath(genpath('/auto/k2/share/strflab'));  % remove current STRF lab path
addpath(genpath('/auto/k2/share/strflabGO'));  % add new STRF lab path

stimEst=x(1:900,:);
respEst=y(1:900);
stimVal=x(901:1000,:);
respVal=y(901:1000);


% Declaring global variables.
%-----------------------------------
global globDat;  % Must declare the global variable globDat in all functions that will access stim and resp.
strfData(stimEst,respEst);


%Initialize strf
%--------------------------------------
strf=linInit(100,[0:4]);

options = trnDirectFit3;
options.time_lag_max = 4;
options.cacheDir = '/auto/fdata/pgill/first_cache_dir';
options.Tol_val = 0;
strf = 1;
datIdx = 2;
[strfTrained,options]=strfOpt(strf,datIdx,options);


% if ~exist('x','var')
%     load fakedata/sampledata1.mat;
% end
% addpath('./lower_train');
% options.cacheDir = '/auto/fdata/pgill/first_cache_dir';% 'null'; %;
% options.Tol_val = .001;
% if exist('assign','var')
%     options.assign = assign;
% end
% 
% if exist('groups','var')
%     options.assign = groups;
% end
% options.time_lag_min = 0;
% options.time_lag_max = 2;
% 
% if exist('hashes','var')
%     options.hashes = hashes;
% end
% stim = x';
% if exist('t','var')
%     resp = t';
% else
%     resp = y';
% end
% 
% 
% [strf, options] = trnDirectFit(stim, resp, options);