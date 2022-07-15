function BCavg_180(verse_factor, display, spokesdefpath, use_fixed_subpulse, data, gDown, opt)
%function BCavg_180(verse_factor, display, spokesdefpath, use_fixed_subpulse)


%% INITIALIZATIONS

close all


%% OPTIONS

opt.do_mls=0;
opt.refoc_pulse=1;
opt.tfa = 180;
opt.bc_pulse = 1;

avg_target = 1;   % for BC pulse, make 1 if you want the average to be the target, 0 if you want the center of the image to be target: scaling of the pulse


%% RUN

%BCpulse(spokesdefpath,display,refoc_pulse,target_phase,verse_factor,b0_zero,avg_target,data,tfa, use_fixed_subpulse, gDown, freqs, freq_weights)
BCpulse(spokesdefpath,opt,avg_target)



%% DISPLAY RESULTS

load all

errorLTA = errorLTA{1};
fprintf('\n ____________ BC ________________\n')
fprintf('\n RMMSE \n')
fprintf('FA: %2.2f %%\n', 100*errorLTA.FArmmse);
fprintf('MT: %2.2f %%\n', 100*errorLTA.MTrmmse);
fprintf('beta2: %2.2f %%\n', 100*errorLTA.betarmmse);
fprintf('\n RMSE \n')
fprintf('MT: %2.2f %%\n', 100*errorLTA.MTrmse);
fprintf('beta2: %2.2f %%\n', 100*errorLTA.betarmse);

