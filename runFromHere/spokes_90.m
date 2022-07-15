function spokes_90(spokesdefpath,opt,lims,i)
%function spokes_90(spokesdefpath,opt)



%% RUN

Spokes_LFA(spokesdefpath,opt,lims,i);



%% PRINTOUT RESULTS

    
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


  
