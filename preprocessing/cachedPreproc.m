function [PS, params] = cashedPreproc(params, stims, fcash);

if ~exist('fcash', 'var')
	 fcash = 1;
end

sdiv1 = ceil(length(stims(:))/929);
sdiv2 = ceil(length(stims(:))/503);
sdiv3 = ceil(length(stims(:))/11);

a.params = params;
a.stimsize = size(stims);
a.stims_fragments1 = stims(1:sdiv1:end);
a.stims_fragments2 = stims(2:sdiv2:end);
a.stims_fragments3 = stims(3:sdiv3:end);


ppdir = params.tmp_dir;
uniquestr = getIDstring(a);
ppfile = [ppdir uniquestr '.mat'];

if exist(ppfile, 'file') & fcash == 1
	disp(['load cashed preprocessed data:' ppfile]);
	load(ppfile);
else
	[PS, params] = feval(params.class, stims, params);
    % [PS, params] = normalizeChannels(PS, params);

    if fcash
        disp(['save preprocessed data:' ppfile]);
        save(ppfile, 'PS', 'params', '-V7.3');
    end
end
