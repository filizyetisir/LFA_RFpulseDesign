function spokes_180(spokesdefpath,opt,lims,i)
%function spokes_180(spokesdefpath,opt)

close all

%% RUN

Spokes_LFA(spokesdefpath,opt,lims,i);



%% DISPLAY RESULTS

load all

N = 1;

errorSTA = errorSTA{N};
fprintf('\n_________RMMSE________________\n')
fprintf('\n STA \n')
fprintf('FA: %2.2f %%\n', 100*errorSTA.FArmmse);
fprintf('MT: %2.2f %%\n', 100*errorSTA.MTrmmse);
fprintf('beta2: %2.2f %%\n', 100*errorSTA.betarmmse);
errorLTA = errorLTA{N};
fprintf('\n  LFA \n')
fprintf('FA  %2.2f %%\n', 100*errorLTA.FArmmse);
fprintf('MT  %2.2f %%\n', 100*errorLTA.MTrmmse);
fprintf('beta2  %2.2f %%\n', 100*errorLTA.betarmmse);
fprintf('\n_________RMSE________________\n')
fprintf('\n STA \n')
fprintf('MT: %2.2f %%\n', 100*errorSTA.MTrmse);
fprintf('beta2: %2.2f %%\n', 100*errorSTA.betarmse);
fprintf('\n  LFA \n')
fprintf('MT  %2.2f %%\n', 100*errorLTA.MTrmse);
fprintf('beta2  %2.2f %%\n\n\n', 100*errorLTA.betarmse);

 

%% SAVE REFOCUSING PHASE for 90


if(opt.create_lcurve == 0)
for i = 1:opt.nfreqs
tmp = angle(-pulse.LFA.BETA2{i});
tmp2 = unwrap_2d_phase_image(tmp, abs(-pulse.LFA.BETA2{i}), 1, 1);
phaseTarget(:,:,i) = (tmp2)/2;
end
% figure; imagesc(phaseTarget); axis image off; colormap(hot); colorbar; caxis([-pi pi])
save phaseTargetFrom180 phaseTarget

end

