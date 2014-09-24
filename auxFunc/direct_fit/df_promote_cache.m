function [outfile,found] = df_promote_cache(the_checksum,cached_dir);
%  This function updates the stats on the usage of a cahsed file.

outfile = fullfile(cached_dir,'metadata',[the_checksum '.mat']);
last_used = now;  %last_used is 1 + the number of days (plus the fraction of the day gone by now) since 01-Jan-0000
%Note: this format is NOT y10k compatible.  Please do
%not sue my great^250-grandchildren for this
%oversight (esp since it's likely you'll be one of
%them...)

if exist(outfile,'file')
    try  %  If another machine tries to update this metadata at the same time (causing an error), who cares?  The correct date is still set.
        save(outfile,'last_used','-APPEND');
    catch
        disp('Had trouble updating a "last_used" variable in the cache metadata.  Is another process trying to update it too?');
        disp('Unless this message appears all the time, it''s not a problem.')
    end
found = 1;
else
    dirout = dir(fullfile(cached_dir,[the_checksum '.mat']));
    space_needed = dirout.bytes / 2^30;  %Space_needed in GB throughout.
    time_saved = .0000012345;
    try  %  If another machine tries to update this metadata at the same time (causing an error), who cares?  The correct date is still set.
    
    save(outfile,'last_used','space_needed','time_saved');
        catch
        disp('Had trouble creating a file in the cache metadata.  Is another process trying to update it too?');
        disp('Unless this message appears all the time, it''s not a problem.')

    end

    found = 0;
end




%  Previous version:
%
% loaded = load(fullfile(cached_dir,'cache_useage_stats.mat'));
% cache_useage_stats = loaded.cache_useage_stats;
% total_space_used = loaded.total_space_used;
% done = 0;
% for jj = 1:length(cache_useage_stats)
%     if ~done
%         if strcmp(the_checksum,cache_useage_stats(jj).checksum)
%             cache_useage_stats(jj).times_used = cache_useage_stats(jj).times_used + 1;
%             done = 1;
%         end
%     end
% end
%
%
% save(fullfile(cached_dir,'cache_useage_stats.mat'),'cache_useage_stats','total_space_used');
