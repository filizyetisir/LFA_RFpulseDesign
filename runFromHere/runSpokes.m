
close all
clear all

mkdir output
path = pwd;
path = path(1:(end-11));
addpath([path, 'Spokes_LFA_v12_filiz_v2']);
addpath([path,'mintverse_v10_filiz']);


%% INITIALIZATION

% Pulse design Options
opt.display = 'on';  % show figures
opt.opt_grad = 0; % spoke optimization
opt.freqs = [0];  % frequency offsets to design the pulse for, for B0 robustness
opt.freq_weights = [1];  % the weights of the frequencies

% For 2 spokes pulses
freqs_2sp = 0;%[0,50,-50]; % Assign to 2 spokes below
freq_weights_2sp = 1;%[10,1,1];
%

% Data
simData = 1; % use simulation data
opt.b0_zero = 0; % Make B0 zero for test purposes
opt.data = 1; % whole 1, inner 2, special 3
opt.downsample = 1; % downsample B1 B0 maps to speed up pulse design
opt.dsr = 2; % downsampling rate, integer

% For ST=3mm, TBW=2.5, gdown=20, subpulse is 1.5 ms 
% and total 2.5 ms for 1spoke, 4.8 ms for 2 spokes
gDown =20;


% For BC center pulses
opt.fil1 = 0;
opt.fil2 = 5;
if(opt.downsample == 1)
    opt.fil1 = round(opt.fil1/opt.dsr);
    opt.fil2 = round(opt.fil2/opt.dsr);
end
%

% Pulse shape and duration
opt.verse_factor = 1; % to overcome peak power limitations and shorten the pulse
opt.gsym = 1; % Make the gradient trajectory completely symmetrical
opt.use_fixed_subpulse = 1; % use fixed subpulse shape instead of calculating it every time
opt.gDown = gDown; % control pulse duration 

% GPU or CPU parallelization to speed up, make only one 1
opt.use_gpu = 1;
opt.multiCPU = 0;
% Jacobian, Hessian
opt.useJ = 1;
opt.useH = 1;
% SAR constraints
opt.sarcons = 1;


% Pulse type, make only one of them 1 for lcurve
bcavg = 0;
bccenter = 0;
onespLS = 0; % For this simulation data 1spLS does not work well, probably due to the unrealistic phase of the CP mode?
onespMLS = 0;
twospMLS = 1;

% CP mode initialization? To avoid the hole solution in 1 spoke case
use_bc_pha = 0;
if(onespMLS == 1 && simData == 1)
    opt.use_bc_init_STA = 0;
    opt.use_bc_init_LFA = 1;
else
    opt.use_bc_init_STA = 0;
    opt.use_bc_init_LFA = 0;
end

% Other
if(use_bc_pha == 1 && onespMLS == 1)
    load bc_pha
    tp_180 = bc_pha;
else
    tp_180 = zeros(250,140);
    % Downsample target phase
    if(opt.downsample == 1)
        tp_180 = tp_180(1:opt.dsr:end, 1:opt.dsr:end);
    end
end


% L curves
opt.create_lcurve = 0;
if(opt.create_lcurve == 1)
    resultsDir = 'lcurves';
else
    resultsDir = 'test';
end

sarlimited = 0;
rflimited = 1;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
duty_cycle = 0.15; % 15%

if(onespLS == 1)
    if(sarlimited == 1)
        all_lsar =10/duty_cycle; % 10 W/kg for peak local SAR limit
        lcurvelims = all_lsar;
    elseif(rflimited == 1)
        all_rf = 300;
        lcurvelims = all_rf;
    end
elseif(onespMLS == 1)
    if(sarlimited == 1)
        if(opt.create_lcurve == 1)
        all_lsar = [25, 27, 28.5, 30, 31.5, 33, 35, 37, 39, 42, 44, 46, 48, 50, 53, 56, 59, 64, 67];
        else
            all_lsar =10/duty_cycle; % 10 W/kg for peak local SAR limit
        end
        lcurvelims = all_lsar;
    elseif(rflimited == 1)
        if(opt.create_lcurve == 1)
            all_rf = [160, 165, 170, 175, 180, 186, 192, 200, 210, 220, 230, 240, 246, 252, 260, 280, 300];
        else
            all_rf = 300;
        end
        lcurvelims = all_rf;
    end
elseif(twospMLS == 1)
    if(sarlimited == 1)
        if(opt.create_lcurve == 1)
            all_lsar = [6.6 7, 7.6, 8.2, 8.8, 9.4, 10, 11, 12, 13, 14, 15, 17, 20, 25, 30, 35, 40, 45, 50,55,  60, 67];
        else
        	all_lsar = 10/duty_cycle; % 10 W/kg for peak local SAR limit
        end
        lcurvelims = all_lsar;
    elseif(rflimited == 1)
        if(opt.create_lcurve == 1)
            all_rf = [83, 86, 89, 92, 95, 98, 101, 105, 110, 115, 120, 125, 130, 140, 150, 160, 180, 200, 220, 240, 270, 300];
        else
            all_rf = 300;
        end
        lcurvelims = all_rf;
    end
elseif(bcavg == 1 || bccenter == 1)
    % BC
    all_lsar = 10/duty_cycle; % 10 W/kg for peak local SAR limit
    lcurvelims = all_lsar;
end


for i =1:length(lcurvelims)
    
    if(sarlimited == 1)
        lims.lsarlim = all_lsar(i);
    else
        lims.lsarlim = 10/duty_cycle;
    end
    if(rflimited == 1)
        lims.rf_max = all_rf(i);
    else
        lims.rf_max = 300;
    end
    lims.gsarlim = 3.2/duty_cycle; % 3.2 W/kg for head average SAR
    lims.avpowlim = 10/duty_cycle; % 10 W per channel average power limit
    
    
%% DESIGN 1 SPOKE LS 180 1 SPOKE LS 90

if(onespLS == 1)
       
    % 180
    optSpokes = opt;
    % Use fixed sub pulse
    if(optSpokes.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    optSpokes.do_mls = 0;
    optSpokes.refoc_pulse = 1;
    optSpokes.tphase = tp_180;
    spokesdefpath = 'spokes_def_180.txt';
    spokes_180(spokesdefpath,optSpokes,lims,i)
    eval(sprintf('mkdir ./Results/%s/1spLS_1spLS/180',resultsDir))
    load all.mat
    if(opt.create_lcurve == 0)
         save(sprintf('./Results/%s/1spLS_1spLS/180/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    else
        save(sprintf('./Results/%s/1spLS_1spLS/180/all%d', resultsDir, i), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA', 'errorSTA', 'MTtimely');
    end
    spcoords = spokes.spcoords;
    
    lcurve(i,1) = sar.lsarmax;
    lcurve(i,2) = sar.gsarmax;
    lcurve(i,3) = sar.ppavmax;
    lcurve(i,4) = sar.ppmax;
    
    lcurve(i,6) = sar.LTA.maxLocSAR_VOP;
    lcurve(i,7) = sar.LTA.globSAR_VOP;
    lcurve(i,8) = max(pulse.LFA.avgp);
    lcurve(i,9) = sum(pulse.LFA.avgp);
    lcurve(i,10) = max(pulse.LFA.maxp);
    
    lcurve(i,12) = errorLTA{1}.betarmmse;
    
    lcurve(i,14) = pulse.STA.comptime;
    lcurve(i,15) = pulse.LFA.timeAll(end);
    
    if(opt.create_lcurve == 0)
        % 90
        optSpokes = opt;
        % Use fixed sub pulse
        if(optSpokes.use_fixed_subpulse)
            load verse_subpulse_1sp.mat
            save('verse_subpulse.mat','rfv','gv')
        end
        %
        optSpokes.do_mls = 0;
        optSpokes.refoc_pulse = 0;
        load phaseTargetFrom180.mat
        optSpokes.tphase= phaseTarget;
        spokesdefpath = 'spokes_def_90.txt';
        spokes_90(spokesdefpath,optSpokes, lims, i)
        eval(sprintf('mkdir ./Results/%s/1spLS_1spLS/90',resultsDir))
        load all.mat
        save(sprintf('./Results/%s/1spLS_1spLS/90/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    end
end


%% DESIGN 1 SPOKE MLS 180 1 SPOKE LS 90

if(onespMLS == 1)
    
    %resultsDir = 'lcurves';
    
    % 180
    optSpokes = opt;
    % Use fixed sub pulse
    if(optSpokes.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    optSpokes.do_mls = 1;
    optSpokes.refoc_pulse = 1;
    optSpokes.tphase = tp_180;
    spokesdefpath = 'spokes_def_180.txt';
    spokes_180(spokesdefpath,optSpokes,lims,i)
    eval(sprintf('mkdir ./Results/%s/1spLS_1spMLS/180',resultsDir))
    load all.mat
    if(opt.create_lcurve == 0)
        save(sprintf('./Results/%s/1spLS_1spMLS/180/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    else
        save(sprintf('./Results/%s/1spLS_1spMLS/180/all%d', resultsDir, i), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA', 'errorSTA', 'MTtimely');
    end
    spcoords = spokes.spcoords;
    
%     xrf = pulse.LFA.xrf;
%     save xrf xrf
    
    lcurve(i,1) = sar.lsarmax;
    lcurve(i,2) = sar.gsarmax;
    lcurve(i,3) = sar.ppavmax;
    lcurve(i,4) = sar.ppmax;
    
    lcurve(i,6) = sar.LTA.maxLocSAR_VOP;
    lcurve(i,7) = sar.LTA.globSAR_VOP;
    lcurve(i,8) = max(pulse.LFA.avgp);
    lcurve(i,9) = sum(pulse.LFA.avgp);
    lcurve(i,10) = max(pulse.LFA.maxp);
    
    lcurve(i,12) = errorLTA{1}.betarmmse;
    
    lcurve(i,14) = pulse.STA.comptime;
    lcurve(i,15) = pulse.LFA.timeAll(end);
    
    if(opt.create_lcurve == 0)
    % 90
    optSpokes = opt;
    % Use fixed sub pulse
    if(optSpokes.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    %
    optSpokes.do_mls = 0;
    optSpokes.refoc_pulse = 0;
    load phaseTargetFrom180.mat
    optSpokes.tphase= phaseTarget;
    spokesdefpath = 'spokes_def_90.txt';
    spokes_90(spokesdefpath,optSpokes, lims, i)
    eval(sprintf('mkdir ./Results/%s/1spLS_1spMLS/90',resultsDir))
    load all.mat
    save(sprintf('./Results/%s/1spLS_1spMLS/90/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    end
end


%% DESIGN 2 SPOKES MLS 180 2 SPOKES LS 90

if(twospMLS == 1)
    
    %resultsDir = 'lcurves';
   
    % 180
    optSpokes = opt;
    % Change sub pulse duration!
    %optSpokes.gDown = 35;
    %
    % Change frequencies 
    optSpokes.freqs = freqs_2sp;%[0 50 -50];
    optSpokes.freq_weights = freq_weights_2sp;%[10 1 1];
    %
    % Use fixed sub pulse
    if(optSpokes.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    %
    optSpokes.do_mls = 1;
    optSpokes.refoc_pulse = 1;
    optSpokes.tphase = repmat(tp_180, [1,1,length(optSpokes.freqs)]);
    spokesdefpath = 'spokes_def_180_2sp.txt';
    spokes_180(spokesdefpath,optSpokes,lims,i)
    eval(sprintf('mkdir ./Results/%s/2spLS_2spMLS/180',resultsDir))
    load all.mat
    if(opt.create_lcurve == 0)
        save(sprintf('./Results/%s/2spLS_2spMLS/180/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    else
        save(sprintf('./Results/%s/2spLS_2spMLS/180/all%d', resultsDir, i), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA','errorSTA', 'MTtimely');
    end
    spcoords = spokes.spcoords;
    
    %xrf = pulse.LFA.xrf;
    %save xrf xrf
    
    lcurve(i,1) = sar.lsarmax;
    lcurve(i,2) = sar.gsarmax;
    lcurve(i,3) = sar.ppavmax;
    lcurve(i,4) = sar.ppmax;
    
    lcurve(i,6) = sar.LTA.maxLocSAR_VOP;
    lcurve(i,7) = sar.LTA.globSAR_VOP;
    lcurve(i,8) = max(pulse.LFA.avgp);
    lcurve(i,9) = sum(pulse.LFA.avgp);
    lcurve(i,10) = max(pulse.LFA.maxp);
    
    lcurve(i,12) = errorLTA{1}.betarmmse;
    
    lcurve(i,14) = pulse.STA.comptime;
    lcurve(i,15) = pulse.LFA.timeAll(end);
    
    
    if(opt.create_lcurve == 0)
        % 90
        optSpokes = opt;
        % Change sub pulse duration!
        %optSpokes.gDown = 35;
        %
        % Change frequencies
        optSpokes.freqs = freqs_2sp;%[0 50 -50];
        optSpokes.freq_weights = freq_weights_2sp;%[10 1 1];
        %
        % Use fixed sub pulse
        if(optSpokes.use_fixed_subpulse)
            load verse_subpulse_1sp.mat
            save('verse_subpulse.mat','rfv','gv')
        end
        %
        optSpokes.do_mls = 0;
        optSpokes.refoc_pulse = 0;
        load phaseTargetFrom180.mat;
        optSpokes.tphase= phaseTarget;
        spokesdefpath = 'spokes_def_90_2sp.txt';
        spokes_90(spokesdefpath,optSpokes,lims,i)
        eval(sprintf('mkdir ./Results/%s/2spLS_2spMLS/90',resultsDir))
        load all.mat
        save(sprintf('./Results/%s/2spLS_2spMLS/90/all_J_%d_H_%d_G_%d_MC_%d', resultsDir, optSpokes.useJ, optSpokes.useH, optSpokes.use_gpu, optSpokes.multiCPU), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorSTA', 'errorLTA', 'MTtimely');
    end
    

end



%% DESIGN BC PULSES 


if(bcavg == 1)
    
    optBCavg = opt;
    optBCavg.do_mls = 0;
    optBCavg.refoc_pulse = 1;
    optBCavg.tphase = tp_180;
    optBCavg.bc_pulse = 1;
    optBCavg.avg_target = 1;
    % Fixed verse subpulse
    if(optBCavg.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    spokesdefpath = 'spokes_def_180.txt';
    BCpulse(spokesdefpath, optBCavg)
    eval(sprintf('mkdir ./Results/%s/BCavg/180',resultsDir))
    load all.mat
    save(sprintf('./Results/%s/BCavg/180/all', resultsDir), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA');
    
    lcurve(i,1) = sar.lsarmax;
    lcurve(i,2) = sar.gsarmax;
    lcurve(i,3) = sar.ppavmax;
    lcurve(i,4) = sar.ppmax;
    
    lcurve(i,6) = sar.maxLocSAR_VOP;
    lcurve(i,7) = sar.globSAR_VOP;
    lcurve(i,8) = max(pulse.LFA.avgp);
    lcurve(i,9) = sum(pulse.LFA.avgp);
    lcurve(i,10) = max(pulse.LFA.maxp);
    
    lcurve(i,12) = errorLTA.betarmmse;
    
    % 90
    optBCavg = opt;
    optBCavg.do_mls = 0;
    optBCavg.refoc_pulse = 0;
    optBCavg.tphase = tp_180;
    optBCavg.bc_pulse = 1;
    optBCavg.avg_target = 1;
    % Fixed verse subpulse
    optBCavg.use_fixed_subpulse = 1;
    if(optBCavg.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    spokesdefpath = 'spokes_def_90.txt';
    BCpulse(spokesdefpath, optBCavg)
    eval(sprintf('mkdir ./Results/%s/BCavg/90',resultsDir))
    load all.mat
    save(sprintf('./Results/%s/BCavg/90/all', resultsDir), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA');
end

if(bccenter == 1)
    
    %resultsDir = 'lcurves';
    
    % 180
    optBCcenter = opt;
    optBCcenter.do_mls = 0;
    optBCcenter.refoc_pulse = 1;
    optBCcenter.tphase = tp_180;
    optBCcenter.bc_pulse = 1;
    optBCcenter.avg_target = 0;
    % Fixed verse subpulse
    if(optBCcenter.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    spokesdefpath = 'spokes_def_180.txt';
    %BCcenter_180(optBCcenter.verse_factor, optBCcenter.display, spokesdefpath, optBCcenter.use_fixed_subpulse, optBCcenter.data, optBCcenter.gDown, tp_180)
    BCpulse(spokesdefpath, optBCcenter)
    eval(sprintf('mkdir ./Results/%s/BCcenter/180',resultsDir))
    load all.mat
    save(sprintf('./Results/%s/BCcenter/180/all', resultsDir), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA');
    
    lcurve(i,1) = sar.lsarmax;
    lcurve(i,2) = sar.gsarmax;
    lcurve(i,3) = sar.ppavmax;
    lcurve(i,4) = sar.ppmax;
    
    lcurve(i,6) = sar.maxLocSAR_VOP;
    lcurve(i,7) = sar.globSAR_VOP;
    lcurve(i,8) = max(pulse.LFA.avgp);
    lcurve(i,9) = sum(pulse.LFA.avgp);
    lcurve(i,10) = max(pulse.LFA.maxp);
    
    lcurve(i,12) = errorLTA.betarmmse;
    
    
    % 90
    optBCcenter = opt;
    optBCcenter.do_mls = 0;
    optBCcenter.refoc_pulse = 0;
    optBCcenter.tphase = tp_180;
    optBCcenter.bc_pulse = 1;
    optBCcenter.avg_target = 0;
    % Fixed verse subpulse
    optBCcenter.use_fixed_subpulse = 1;
    if(optBCcenter.use_fixed_subpulse)
        load verse_subpulse_1sp.mat
        save('verse_subpulse.mat','rfv','gv')
    end
    spokesdefpath = 'spokes_def_90.txt';
    %BCcenter_90(optBCcenter.verse_factor, optBCcenter.display, spokesdefpath, optBCcenter.use_fixed_subpulse, optBCcenter.data, optBCcenter.gDown, tp_180)
    BCpulse(spokesdefpath, optBCcenter)
    eval(sprintf('mkdir ./Results/%s/BCcenter/90',resultsDir))
    load all.mat
    save(sprintf('./Results/%s/BCcenter/90/all', resultsDir), 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA');
end

end

save lcurve lcurve
