function out= df_initialize_dir_for_caching(cached_dir);
% cache_useage_stats = struct([]);
% total_space_used = 0;
% save(fullfile(cached_dir,'cache_useage_stats.mat'),'cache_useage_stats','
% total_space_used');
df_rec_make_dir(fullfile(cached_dir,'metadata'));
