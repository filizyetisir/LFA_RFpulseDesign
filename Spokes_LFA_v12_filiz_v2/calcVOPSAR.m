function sar = calcVOPSAR(sar,spokes,system,rf)
%function sar = calcVOPSAR(sar,spokes,system,rf)



%% INITIALIZE

L = length(rf);
xc = rf(1:L/2)+1i*rf((L/2+1):end);
nSAR = size(sar.sarmats,1);
cons = spokes.sumsincsq/spokes.ntimes;
xs = zeros(spokes.nspokes,system.ncoils);
for sp = 1:spokes.nspokes
    xs(sp,:) = xc((sp-1)*system.ncoils+1:sp*system.ncoils);
end




%% CALCULATE SAR


% Local and Global SAR
sarvalues = zeros(nSAR,1);
for n = 1:nSAR
    for sp = 1:spokes.nspokes
        sarvalues(n) = sarvalues(n) + xs(sp,:)*squeeze(sar.sarmats(n,:,:))*xs(sp,:)';
    end
end
sarvalues = sarvalues*cons;
[maxLocalSar, maxLocalSarInd] = max(sarvalues(1:(nSAR-1)));
globSAR = sarvalues(nSAR);


% Local SAR Overestimation
overEst = 0;
for sp = 1:spokes.nspokes
   	overEst = overEst + xs(sp,:)*sar.umax*xs(sp,:)';
end
overEst = overEst*cons;


% Output assignment
sar.maxLocSAR_VOP = maxLocalSar;
sar.maxLocSARInd_VOP = maxLocalSarInd;
sar.globSAR_VOP = globSAR;
sar.lSARoverEst = overEst;




end

