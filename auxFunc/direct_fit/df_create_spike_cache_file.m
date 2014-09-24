function hashes_of_spikes = df_create_spike_cache_file(outputPath,DS);
%  Creates the file with the hashes of the stimuli
hashes_of_spikes = {};
df_rec_make_dir(fullfile(outputPath,'spike_hashes'));
for ii = 1:length(DS)
    dsname = DS{ii}.respfiles;
    [dsdir,dsname,dsext] = fileparts(dsname);
    dsname = [dsname dsext];
    filename = fullfile(outputPath,'spike_hashes',dsname);
    if ~exist(filename,'file')
        this_hash = df_checksum_from_file(DS{ii}.respfiles);
        save(filename,'this_hash');
    else
        load(filename);
    end
    hashes_of_spikes{ii} = this_hash;
end
