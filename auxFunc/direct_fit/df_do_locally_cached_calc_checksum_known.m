function varargout = df_do_locally_cached_calc_checksum_known(local_cache_dir,function_name,the_checksum,varargin)
%  Returns function_name(varargin), using a cache of pre-computed answers when possible.
%  This function is more faster than do_cached_calc but
%  isn't foolproof.  Use do_cached_calc_checksum_known if you know caching
%  is enabled and if the inputs are big in terms of memory/data.  It takes
%  a few hendred ms to make a checksum from 1 MB of data, so if you're
%  going to be re-useing the same data as inputs to different functions,
%  it's better to re-use a *checksum* of the input instead of the input itself
%  when figuring out a unique checksum depending on the data.  (Hence the
%  reason for calling do_cached_calc_checksum_known instead.)

[cached_dir,maxsize] = df_dir_of_caches;
cached_dir = local_cache_dir;
if strcmp(cached_dir,'null')
    disp('Caching is currently disabled.  To learn more about caching, type "help Caching_Readme".');
    %Why does MATLAB choke on "varargout = feval(function_name,varargin)"?
    %With this limitation, I can't even make this hack modular.
    %This next bit is a hack which gets around this limitation.
    if nargout == 0
        to_eval = '';
    else
        to_eval = '[';
        for jj = 1:nargout
            to_eval = [to_eval 'varargout{' num2str(jj) '},'];
        end
        to_eval(end) = ']';
        to_eval = [to_eval ' = '];
    end
    to_eval = [to_eval 'feval(function_name'];
    for jj = 1:(nargin-3)
        to_eval = [to_eval ',varargin{' num2str(jj) '}'];
    end
    to_eval = [to_eval ');'];
    tic;
    eval(to_eval);
    time_taken = toc;
else
    if ~exist(cached_dir,'dir')
        disp(['Message from function "df_do_locally_cached_calc_checksum_known": making cache directory "' cached_dir '".']);
        df_rec_make_dir(cached_dir);
    end
    df_initialize_dir_for_caching(cached_dir);

    the_filename = fullfile(cached_dir,[the_checksum '.mat']);
    need_to_calculate = 1;
    if exist(the_filename,'file')
        try
            loaded = load(the_filename);
            if ~isfield(loaded,'cached_no_args') %Backwards compatability w/ a spelling mistake
                loaded.cached_no_args = loaded.cashed_no_args;
                cached_no_args = loaded.cashed_no_args;
                save(filename,'cached_no_args','-APPEND');
            end
            if loaded.cached_no_args >= nargout
                %disp(['Score one for caching: we don''t have to evaluate the function ' function_name ' since its output has been cached.']);
                for jj = 1:nargout
                    eval(['varargout{jj} = loaded.out' num2str(jj) ';']);
                end
                need_to_calculate = 0;
                [junk,found] = df_promote_cache(the_checksum,cached_dir);
                if ~found
                    disp(['It''s funny: I just used the data for the calculation ' the_checksum char(10) ...
                        ', but I don''t know how much time it saved.  Assigning flag value of .0000012345 to time_saved.']);
                    % I would use -999 as the flag value, but it really should be positive but small.
                end


            end
        catch
            disp(['File ' the_filename ' was corrupt.  Deleting...']);
            try
                delete(the_filename)
            catch
                disp(lasterr)
            end
        end

    end

    if need_to_calculate == 1
        %Why does MATLAB choke on "varargout = feval(function_name,varargin)"?
        %With this limitation, I can't even make this hack modular.
        %This next bit is a hack which gets around this limitation.
        if nargout == 0
            to_eval = '';
        else
            to_eval = '[';
            for jj = 1:nargout
                to_eval = [to_eval 'varargout{' num2str(jj) '},'];
            end
            to_eval(end) = ']';
            to_eval = [to_eval ' = '];
        end
        to_eval = [to_eval 'feval(function_name'];
        for jj = 1:(nargin-3)
            to_eval = [to_eval ',varargin{' num2str(jj) '}'];
        end
        to_eval = [to_eval ');'];
        tic;
        eval(to_eval);
        time_taken = toc;
        %Hack finished.
        %Now a quick benchmark to keep time_taken meaningful:
        %tic; A = inv(rand(500)); t_norm = toc;%  To normalize for processor speed
        %tic; A = inv(rand(500)); t_norm = toc;%  since loading the inv and rand functions might take time
        t_norm = 1;  % Do we really need to benchmark?  Maybe it's just useless overhead.
        time_saved = time_taken/t_norm;
        %%% Now let's see if there's enough space to save the file  We're going
        %%% to use the following procedure to figure out if we should save the
        %%% file.  If cached_dir is > 10% empty and there's space for the new
        %%% file, it gets saved.  If not, we'll clean up the cached_dir so it's
        %%% at least 20% empty (taking the most useless files away first), then
        %%% save the file if there's room for it.
        
        whos_varargout = whos('varargout');
        try
            space_needed = whos_varargout.bytes / 2^30;  %  For size in GB - won't be exactly file size, but close enough for now
        catch
            space_needed = 1e-9;
        end
        total_space_used = df_get_total_space_used(cached_dir);
        if (total_space_used < (.9*maxsize)) & ((space_needed + total_space_used) < maxsize)
            willsave = 1;
        else
            if total_space_used > (.8*maxsize)
                desired_cache_size = .8*maxsize;
                try
                    do_cache_cleanup(cached_dir,desired_cache_size);

                catch
                    disp(lasterr)
                    pause(.1+rand(1))
                    try
                        do_cache_cleanup(cached_dir,desired_cache_size);
                    catch
                        disp(lasterr)

                        disp('Trying to clean up the cache, but it looks like another computer is trying too.  I give up.')
                    end
                end

            end
            if space_needed < (.2*maxsize)
                willsave = 1;
            else
                willsave = 0;
            end
        end
        if willsave
            for jj = 1:nargout
                eval(['out' num2str(jj) ' = varargout{' num2str(jj) '};']);
            end
            cached_no_args = nargout;
            cashed_no_args = cached_no_args;
            try
                if ~exist(the_filename,'file')
                    disp(['Saving the results from evaluating function ' function_name ' in a local directory for future use.'])
                    save(the_filename,'out*','cached_no_args','cashed_no_args');  %Backwards compatability w/ a spelling mistake
                    df_add_record_to_cache(cached_dir,the_checksum,time_saved,function_name,the_filename);
                else
                    t_loaded = load(the_filename);
                    if isfield(t_loaded,['out' num2str(cached_no_args)])
                        disp(['I was going to save the results from evaluating function ' function_name ', but another process has beatme to it.']);
                    else
                        disp(['Although function ' function_name ' has been evaluated before with exactly the same inputs, now more outputs have been requested, so I''m caching them all.']);
                        save(the_filename,'out*','cached_no_args','cashed_no_args');  %Backwards compatability w/ a spelling mistake
                        df_add_record_to_cache(cached_dir,the_checksum,time_saved,function_name,the_filename);
                    end
                end
            catch
                disp(['There was a problem saving results from the function ' function_name ' to the cache.  Skipping...']);
                disp(lasterr);
            end

        else
            disp(['It is not worth saving the output of ' function_name ' this time: it would take up too much space.']);

        end

    end
end
