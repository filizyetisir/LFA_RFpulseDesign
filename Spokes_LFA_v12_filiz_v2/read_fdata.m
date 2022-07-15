function y=read_fdata(file,n)
% function y=read_fdata(file,n). Arguments:
% file: input file name
% n: # of data point to read (-1 to read the whole file)

fp=fopen(file,'rb');
if(n>0)
    y=fread(fp,n,'float32');
else
    y=fread(fp,'float32');
end
fclose(fp);
