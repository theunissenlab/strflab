function df_add_record_to_cache(cached_dir,the_checksum,time_saved,function_name,data_filename);
%  Creates a file to keep track of how useful a calculation was.
%  Also puts the same inofrmation in the file data_filename, which is where
%  the (long) output of the calculation is stored as well.

outfile = df_promote_cache(the_checksum,cached_dir);  %If it doesn't exist, creates a file to keep track of how useful this calculation is.
save(outfile,'time_saved','function_name','-APPEND');
load(outfile);
save(data_filename,'time_saved','function_name','space_needed','-APPEND');

%  Old version:
%function active = df_add_record_to_cache(cached_dir,maxsize,the_checksum,time_saved,space_needed);
%
% %  Adds a record to the cache, or if it's already there, increments it.
% loaded = load(fullfile(cached_dir,'cache_useage_stats.mat'));
% cache_useage_stats = loaded.cache_useage_stats;
% total_space_used = loaded.total_space_used;
% 
% % See if the checksum already exists
% 
% done = 0;
% active = -999;  %  Will be the number of cache_useage_stats corresponding to the file we're looking at.
% for jj = 1:length(cache_useage_stats)
%     if ~done
%         if strcmp(the_checksum,cache_useage_stats(jj).checksum)
%             cache_useage_stats(jj).times_used = cache_useage_stats(jj).times_used + 1;
%             cache_useage_stats(jj).time_saved = time_saved;
%             cache_useage_stats(jj).space_needed = space_needed;
%             done = 1;
%             active = jj;
%         end
%         cache_useage_stats(jj).bargain = sqrt(cache_useage_stats(jj).times_used) * cache_useage_stats(jj).time_saved / cache_useage_stats(jj).space_needed;
%     end
% end
% 
% %  Create entry if needed.
% 
% if ~done
%     disp(['Creating entry in list for checksum ' the_checksum '.']);
%     jj = length(cache_useage_stats) + 1;
%     cache_useage_stats(jj).checksum = the_checksum;
%     cache_useage_stats(jj).times_used = 1;
%     cache_useage_stats(jj).time_saved = time_saved;
%     cache_useage_stats(jj).space_needed = space_needed;
%     cache_useage_stats(jj).bargain = sqrt(cache_useage_stats(jj).times_used) * cache_useage_stats(jj).time_saved / cache_useage_stats(jj).space_needed;
%     active = jj;
% end
% save(fullfile(cached_dir,'cache_useage_stats.mat'),'cache_useage_stats','total_space_used');
