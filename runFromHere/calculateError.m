function error = calculateError(MT, FA, BETA2, adj, opt, i)
%function error = calculateError(MT, FA, BETA2, adj, opt)

ind = find(adj.roi);


target = adj.roi.*opt.tfa;
error.FArmmse = norm(abs(FA(ind))-abs(target(ind)))/norm(abs(target(ind)));

%for i = 1:opt.nfreqs
    [Nx Ny] = size(adj.roiWhole);
    target = adj.targetImg_LFA((i-1)*Nx+1:i*Nx,:);
    %target = adj.roi.*sin(opt.tfa/180*pi).*exp(1i*opt.tphase(:,:,i));
    error.MTrmse = norm(MT(ind)-target(ind))/norm(target(ind));
    error.MTrmmse = norm(abs(MT(ind))-abs(target(ind)))/norm(abs(target(ind)));
%end

% for i = 1:opt.nfreqs
    [Nx Ny] = size(adj.roiWhole);
    target = adj.targetImg_LFA((i-1)*Nx+1:i*Nx,:);
    error.betarmse = norm(BETA2(ind)-target(ind))/norm(target(ind));
    error.betarmmse = norm(abs(BETA2(ind))-abs(target(ind)))/norm(abs(target(ind)));
% end

end

