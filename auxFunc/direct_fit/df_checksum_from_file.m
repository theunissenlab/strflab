function out = df_checksum_from_file(filename);
%  Produces a checksum from a file that's already there.

if exist('/usr/bin/openssl','file')
    [j1,j2] = system(['tail -c +120 ' filename ' | /usr/bin/openssl md5']);%  Remove the date stamp in the file; we don't want our hashes to depend on that!
    out = upper(j2(1:32));
else


    fid = fopen(filename);
    input = fread(fid,inf,'uint8');
    fclose(fid);
    input = uint8(input(120:end));   %  Remove the date stamp in the file; we don't want our hashes to depend on that!
    filter = [8    28   -84     1    69   114   -45   -49   107   -55    92   118   -55 -53  -102   -94];%  Note: there was a problem with the old vector: [ones(16,1) ; -1*ones(16,1)];

    the_conv = conv(filter,double(input));

    out =[];
    for jj = 1:16
        toadd = mod(sum(the_conv(jj:16:end)),256);
        out = [out dec2hex(toadd,2)];
    end
end