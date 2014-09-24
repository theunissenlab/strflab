function hashes_of_stims = df_create_stim_cache_file(outputPath,DS)
%  Creates the file with the hashes of the stimuli
hashes_of_stims = {};
df_rec_make_dir(fullfile(outputPath,'stim_hashes'));
running_flag = 1;
using_bar = 0;
for ii = 1:length(DS)
    dsname = DS{ii}.stimfiles;
    [dsdir,dsname,dsext] = fileparts(dsname);
    dsname = [dsname dsext];
    filename = fullfile(outputPath,'stim_hashes',dsname);
    if ~exist(filename,'file')
        if (running_flag == 1) & using_bar == 0
            using_bar = 1;
        end

        this_hash = df_checksum_from_file(DS{ii}.stimfiles);
        save(filename,'this_hash');
    else
        load(filename);
    end
    hashes_of_stims{ii} = this_hash;
end
