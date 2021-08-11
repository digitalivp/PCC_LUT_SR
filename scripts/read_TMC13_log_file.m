function bpp = read_TMC13_log_file(log_filename)

log_file = fileread(log_filename);
key_V = 'positions bitstream size ';
key_C = 'colors bitstream size ';
key_T = 'Total bitstream size ';

i_V = strfind(log_file,key_V);
i_C = strfind(log_file,key_C);
i_T = strfind(log_file,key_T);
bytes_V = 0;
bytes_C = 0;
for k=1:length(i_V)
  bytes_V = bytes_V + sscanf(log_file(i_V(k)+length(key_V):end),'%g',1);
end
for k=1:length(i_C)
  bytes_C = bytes_C + sscanf(log_file(i_C(k)+length(key_C):end),'%g',1);
end
bytes_T = sscanf(log_file(i_T(1)+length(key_T):end),'%g',1);

bpp   = struct;
bpp.bitstream_V = 8 * bytes_V;
bpp.bitstream_C = 8 * bytes_C;
bpp.bitstream_T = 8 * bytes_T;

end
