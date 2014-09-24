function out = df_rec_make_dir(dirname)
% Creates a directory even if its subdirectories do not exist.

posslash = findstr(dirname,filesep);
for jj = 2:length(posslash)
    to_check = dirname((posslash(jj-1)+1):posslash(jj));
    parent = dirname(1:posslash(jj-1));
    if ~exist(fullfile(parent,to_check),'dir')
        mkdir(parent,to_check);
    end
end
if ~exist(dirname,'dir')
    to_check = dirname((posslash(end)+1):end);
    parent = dirname(1:posslash(end));
    if ~exist(fullfile(parent,to_check),'dir')
        mkdir(parent,to_check);
    end
end
out = 1;
