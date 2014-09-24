%% Gets the directory in which a function resides
function dname = get_function_dir(funcName)

    dname = '';
    fpath = which(funcName);
    if ~isempty(fpath)
        dname = fileparts(fpath);
    end

    