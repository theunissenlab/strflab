function total_space_used = df_get_total_space_used(cached_dir)
dirout = dir(cached_dir);
total_space_used = sum([dirout(:).bytes])/2^30;
