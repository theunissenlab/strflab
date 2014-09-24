function [cached_dir,maxsize] = df_dir_of_caches
%  Here, change 'null' into the full path directory where you want cached
%  results to be saved.  Remember to include the full path, starting with
%  'C:\' for Windows and '/' for other operating systems.

cached_dir = 'null';
%cached_dir = 'C:\STRFPAK_data\first_cache_dir';
%cached_dir = '/Applications/MATLAB7/work/Theunissen/half_hash/first_cache_dir';
%cached_dir = '/auto/fdata/pgill/first_cache_dir';
maxsize = 10;  % target maximum size of the files to put in this directory, in GB.  Give until it hurts.
