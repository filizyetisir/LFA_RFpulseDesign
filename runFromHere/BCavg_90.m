function BCavg_90(verse_factor, display, spokesdefpath, use_fixed_subpulse, data, gDown, tp)
%function BCavg_90(verse_factor, display, spokesdefpath,use_fixed_subpulse)


%% INITIALIZATIONS

%close all


%% OPTIONS


do_mls=0;

refoc_pulse=0;

avg_target = 1;   % for BC pulse, make 1 if you want the average to be the target, 0 if you want the center of the image to be target: scaling of the pulse

target_phase = tp;

opt_grad = 0;


b0_zero = 0; % make it 1 if you want 0 b0_map, for test purposes

use_gpu = 1;

freqs=0;
freq_weights=1;

%data = 2; % whole if 1, inner if 2, special if 3 (mask of the prescan data)

tfa = 90;

%% RUN

BCpulse(spokesdefpath,display,refoc_pulse,target_phase,verse_factor,b0_zero,avg_target,data,tfa, use_fixed_subpulse, gDown, freqs, freq_weights)



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



