

function [c ceq dc dceq]=fn_SP(x,adj,spokes,sar,system)
% function [c ceq dc dceq]=fn_SP(x,adj,spokes,sar,system)

xrf=x( 2*spokes.nspokes+1:end,1 );
[c ceq dc dceq]=fn(xrf,adj,spokes,sar,system);
c=[zeros(2*spokes.nspokes,1) ; c];

dc=[ zeros(2*spokes.nspokes,size(c,1)) ; [zeros(size(dc,1),2*spokes.nspokes) dc] ];




