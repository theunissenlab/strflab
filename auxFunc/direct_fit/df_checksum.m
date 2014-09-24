function out = df_checksum(varargin)
%  Produces a checksum of the input arguments, which can be anything.
%  Note: this is NOT a hash; it does not have crypto strength.  BUT it's
%  still pretty sensitive on its inputs, so it will do for our purposes.

if exist('/usr/bin/openssl','file')
    out = df_checksum_with_md5(varargin);
else

    input = df_concat_for_checksum(varargin);

    filter = [8    28   -84     1    69   114   -45   -49   107   -55    92   118   -55 -53  -102   -94];%  Note: there was a problem with the old vector: [ones(16,1) ; -1*ones(16,1)];

    the_conv = conv(filter,double(input));

    out =[];
    for jj = 1:16
        toadd = mod(sum(the_conv(jj:16:end)),256);
        out = [out dec2hex(toadd,2)];
    end
end