

function [avgp, maxp] = display_optimization_summary(rf,adj,spokes,sar,system)
% function display_optimization_summary(rf,adj,spokes,sar,system)

ncunk=size(rf,1);
nchannels=system.ncoils
xrf=[real(rf);imag(rf)];
% [c ceq dc dceq]=fn(xrf,adj,spokes,sar,system);

% % SAR
% [lsar vopind]=max( c(1:sar.nvop,1)+sar.lsarmax );
% gsar=c(sar.nvop+1,1) + sar.gsarmax;
% disp( sprintf('The estimated local SAR for this set of pulses is %e W/kg (VOP #%d)',lsar,vopind) );
% disp( sprintf('Gocal SAR for this set of pulses is %e W/kg',gsar) );

% power, this is 2 if B1 maps come from experimental set, 8 if from
% simulated data set, but now if in vivo I divide the b1 maps by 2, so it
% can be kept this way, 8
pp=reshape( abs(rf(:,1)).^2/(8.0*50.0),nchannels,spokes.nspokes );
spokesindoff=0;
for i=1:nchannels
    avgp(i) = sum(pp(i,:),2)*(spokes.sumsincsq/spokes.ntimes);
    maxp(i) = max(pp(i,:));
    fprintf('Chn. #%d:\tav. power=%f W\tmax. power=%f W\n',i,sum(pp(i,:),2)*(spokes.sumsincsq/spokes.ntimes),max(pp(i,:)));
end    
fprintf('Total av. power=%f W\t max. av. power / chn=%f W\n', sum(avgp), max(avgp));


