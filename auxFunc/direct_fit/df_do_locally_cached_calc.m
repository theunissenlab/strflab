function varargout = df_do_locally_cached_calc(local_cache_dir,function_name,varargin)
%  Returns function_name(varargin), using a cache of pre-computed answers when possible.
%  This function is more foolproof than do_cached_calc_checksum_known but
%  can be slower.  Use do_cached_calc_checksum_known if you know caching
%  is enabled and if the inputs are big in terms of memory/data.  It takes
%  a few hendred ms to make a checksum from 1 MB of data, so if you're
%  going to be re-useing the same data as inputs to different functions,
%  it's better to re-use a *checksum* of the input instead of the input itself
%  when figuring out a unique checksum depending on the data.  (Hence the
%  reason for calling do_cached_calc_checksum_known instead.)

%  First off, let's see if this computer knows about the cached directory.

if ~exist('df_dir_of_caches.m','file')
    disp('Function "do_cached_calc" cannot find the file "df_dir_of_caches.m", which tells it where to put its saved intermediate results.  Evaluating anyways...');
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
    for jj = 1:(nargin-2)
        to_eval = [to_eval ',varargin{' num2str(jj) '}'];
    end
    to_eval = [to_eval ');'];
    eval(to_eval);
    %Hack finished.
else
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
        for jj = 1:(nargin-2)
            to_eval = [to_eval ',varargin{' num2str(jj) '}'];
        end
        to_eval = [to_eval ');'];
        eval(to_eval);
        %Hack finished.
    else
        if ~exist(cached_dir,'dir')
            disp(['Message from function "do_locally_cached_calc": making cache directory "' cached_dir '".']);
            df_rec_make_dir(cached_dir);
            df_initialize_dir_for_caching(cached_dir);
        end
        if ~exist(fullfile(cached_dir,'metadata'),'dir')
            df_initialize_dir_for_caching(cached_dir);
        end


        %  First step: figure out what the hash of the computation is.

        %  Let's include the actual code going into a function inthe hash itself;
        %  that way if you make any changes to the code, the hash will be
        %  re-computed.  BEWARE: I'm not recursively checking whether the functions
        %  "function_name" calls have also changed since the last hash computation;
        %  that's a risk we'll have to take.

        function_text = df_load_function_text(function_name);

        %b = which('get_strfpak_version');
        if 0%length(b) > 0
            strfpak_version = df_get_strfpak_version;
        else
            strfpak_version = 'unknown';  %Keeps things from crashing if strfpak_version unavailable
        end

        %  Now let's make a hash of the function text and its input, as well as the
        %  strfpak version.
        global the_checksum
        the_checksum = df_checksum(function_text,strfpak_version,varargin);
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
        to_eval = [to_eval 'df_do_locally_cached_calc_checksum_known(local_cache_dir,function_name,the_checksum'];
        for jj = 1:(nargin-2)
            to_eval = [to_eval ',varargin{' num2str(jj) '}'];
        end
        to_eval = [to_eval ');'];
        eval(to_eval);
        %Hack finished.
    end
end
