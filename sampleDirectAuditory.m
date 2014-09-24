if ~exist('x','var')
    load fakedata/sampleDataAuditory.mat;
end
addpath('./lower_train');
% For best results, make this a directory you can write to.
%options.cacheDir = '/auto/fdata/pgill/first_cache_dir';% 'null'; %;
options.cacheDir = '/auto/k5/wafting/cacheDir';% 'null'; %;
options.Tol_val = .001;
if exist('assign','var')
    options.assign = assign;
end

if exist('groups','var')
    options.assign = groups;
end
options.time_lag_min = 2;
options.time_lag_max = 52;

if exist('hashes','var')
    options.hashes = hashes;
end
stim = x;
if exist('t','var')
    resp = t;
else
    resp = y;
end


[strf, options] = trnDirectFit(stim, resp, options);