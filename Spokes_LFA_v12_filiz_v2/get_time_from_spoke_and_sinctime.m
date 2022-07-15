

function time=get_time_from_spoke_and_sinctime(spokes,spokenumber,sinctime)
% function time=get_time_from_spoke_and_sinctime(spokes,spokesnumber,sinctime)

if spokenumber<=0 || spokenumber>spokes.nspokes || sinctime<=0 || sinctime>spokes.ntimesperspoke
    time=-1;
else
    time=spokes.spokes_time_offsets(spokenumber)+sinctime;
end

