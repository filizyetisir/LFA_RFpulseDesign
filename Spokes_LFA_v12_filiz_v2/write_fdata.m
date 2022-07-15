

function write_fdata(data,filepath)
% function write_fdata(data,filepath)

file=fopen(filepath,'wb');
fwrite(file,data,'float32');
fclose(file);