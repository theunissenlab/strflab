function out = df_load_function_text(function_name);
function_path = which(function_name);
fid = fopen(function_path);
if fid > 0
    out = char(fread(fid,inf,'char')');
    fclose(fid);
else
    out = function_name;
end
