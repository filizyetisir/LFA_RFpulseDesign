


function [param1 param2 param3]=tline_to_param(tline,istart,istop)
% function [param1 param2 param3]=tline_to_param(tline,istart,istop)

tmp=tline(istart:istop);

spacechar=sprintf(' ');
tabchar=sprintf('\t');

istart=1; iend=1;
while tmp(iend)~=spacechar && tmp(iend)~=tabchar
    iend=iend+1;
end
param1=str2num(tmp(istart:iend));
while tmp(iend)==spacechar || tmp(iend)==tabchar
    iend=iend+1;
end

istart=iend;
while tmp(iend)~=spacechar && tmp(iend)~=tabchar
    iend=iend+1;
end
param2=str2num(tmp(istart:iend));
while tmp(iend)==spacechar || tmp(iend)==tabchar
    iend=iend+1;
end

istart=iend;
param3=str2num(tmp(istart:end));




