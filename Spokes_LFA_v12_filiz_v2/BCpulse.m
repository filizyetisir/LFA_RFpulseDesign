function BCpulse(filepath,opt)
%function BCpulse(filepath,display,refoc_pulse,tphase,verse_factor,b0_zero,avg_target,data,tfa,use_fixed_subpulse, gDown, freqs, freq_weights)
%function BCpulse(filepath,display,refoc_pulse,tphase,verse_factor,b0_zero,avg_target,data,tfa,use_fixed_subpulse)


%% PREPARE 

% Read parameters
fprintf('Reading pulse parameters ...\n');
[adj, system, opt, spokes, sar]=read_spokes_def_file(filepath,opt);
spokes.SS_grad_down = opt.gDown/100;


% Make b0 zero for test purposes
if(opt.b0_zero == 1)
    adj.b0mapWhole = adj.b0mapWhole*0;
    adj.b0mapInner = adj.b0mapInner*0;
    adj.b0mapSpecial = adj.b0mapSpecial*0;
end

% Other parameters
% opt.refoc_pulse=refoc_pulse;
% opt.tphase=tphase;
% opt.data = data;
% opt.do_mls = 0;
% opt.tfa = tfa;
spokes.verse_factor=opt.verse_factor;
% opt.use_fixed_subpulse = use_fixed_subpulse;
% opt.bc_pulse = 1;
% opt.freqs= freqs;
% opt.freq_weights = freq_weights;
opt.nfreqs = length(opt.freqs);
% Make gradient symmetrically balanced
%opt.gsym=1;

% Check spoke
if(spokes.nspokes > 1)
    error('BC mode pulse can not have more than 1 spoke');
end
if(max(spokes.spcoords)>0)
    error('BC mode pulse can not have a spoke that is not at DC');
end


% Read sar matrices
%fprintf('Reading SAR matrices ...\n');
sar=prepare_sar(sar,spokes,system);


% Prepare optimization
fprintf('Preparing optimization ...\n');
adj=prepare_optimization(opt,spokes,system,adj);
spokes=prepare_pulse(spokes,opt,system);  % prepare pulse on the system raster time (used to export the ini files)


% % Make b0 zero for test purposes
% if(b0_zero == 1)
%     adj.b0map = adj.b0map*0;
%     adj.b0map_ = adj.b0map_*0;
% end


system_opt=system;
spokes_opt=spokes;
spokes_opt.q_gblips=get_Q_gradient_blips(spokes_opt,adj,system_opt);



%% CALCULATE PULSE

[ind1, ind2] = size(adj.roiWhole);
    

%%%% fil1 and fil2 depend on the data
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fil1 = opt.fil1; fil2 = opt.fil2;
% give a 2 cm range
range = round(0.01/adj.dx);
ind1 = round(ind1/2)+fil1; ind2 = round(ind2/2)+fil2;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% BC mode B1+ map
C = system.ncoils;
bcmode = zeros(size(squeeze(adj.b1maps(1,:,:))));
%pha = [0.78; 1.57 ;2.35 ;3.14; 3.93; 4.71 ;5.50; 6.28];
for c = 1:C
    bcmode = bcmode + squeeze(adj.b1maps(c,:,:).*exp(-1j*angle(adj.b1maps(c,ind1,ind2))));
    %bcmode = bcmode + squeeze(adj.b1maps(c,:,:));%*exp(-1j*(pha(c)));
end
figure; 
subplot(1,2,1); imagesc(abs(bcmode)); axis image off; colorbar; colormap(hot); set(gca, 'FontSize', 18);
subplot(1,2,2); imagesc(angle(bcmode)); axis image off; colorbar; colormap(hot); set(gca, 'FontSize', 18);


% Find the scaling factor to achieve the target
gamma = 42.576e6;
if(opt.avg_target == 0) % The center of the FA map will be target FA
    % Small box around the center
    ind1_range = [ind1-range:ind1+range];
    ind2_range = [ind2-range:ind2+range];
    tmp = abs(bcmode(ind1_range,ind2_range));
    tmp = mean(tmp(:));
%     [ind1, ind2] = size(bcmode);
%     ind1 = ind1/2+fil1; ind2 = ind2/2+fil2;
%     while(abs(bcmode(ind1,ind2)) == 0)
%         ind1 = ind1+1;
%         ind2 = ind2+1;
%     end
    scaleFac = opt.tfa/(360*system.deltat*gamma*sum(spokes.sinc)*tmp);
else % The average of the FA map will be target FA
    tmp = (adj.roi).*360*system.deltat*gamma*sum(spokes.sinc).*abs(bcmode);
    avg_FA = sum(tmp(:))/length(find(adj.roi));
    scaleFac = opt.tfa/avg_FA;
end
if(scaleFac > system.umax)
    scaleFac = system.umax;
    error('BC mode RF pulse exceeds maximum voltage limit. Make the pulse duration longer')
end


% Form the BC mode pulse
pha = angle(adj.b1maps(:,ceil(end/2+fil1),ceil(end/2+fil2)));
if(opt.refoc_pulse == 0)
    rf = scaleFac*ones(8,1).*exp(-1j*pha);
else
    rf = scaleFac*ones(8,1).*exp(-1j*pi/2).*exp(-1j*(pha));
end
% for c = 1:8
%     if(opt.refoc_pulse == 0)
%         rf(c,1) = exp(1j*2*pi/8*(c))*scaleFac;
%     else
%         rf(c,1) = exp(1j*(pi/2+2*pi/8*(c)))*scaleFac;
%     end
% end

% Compansate for the slice profile by a scale factor

xrf = [real(rf); imag(rf)];

sar = calcVOPSAR(sar,spokes,system,xrf);


%%  SIMULATE PULSE

rf_full=write_and_display_pulse(rf,opt,spokes,system,'LFA');
if exist('pTXArbitrary.ini','file')
    movefile('pTXArbitrary.ini','pTXArbitrary_LTA.ini');
else
    warning('The LFA pulse file pTXArbitrary.ini does not exist.');
end
if opt.refoc_pulse==1
    M0=[0;1;0];
else
    M0=[0;0;1];
end
nsims=1;  % must be odd to simulate b0_off=0 Hz (Larmor)
b0_offsets=0; % linspace(-100,100,nsims)';
MT=cell(nsims,1);
MZ=cell(nsims,1);
FA=cell(nsims,1);
BETA2=cell(nsims,1);
errorLTA=cell(nsims,1);
for i=1:size(b0_offsets,1)   
    [mx, my, mz, BETA2{i}]=bloch_simulation(rf_full,adj,spokes,system,b0_offsets(i),M0);
    MT{i}=mx + 1j*my;
    MZ{i}=mz;    
    FA{i}=zeros(size(adj.roi));    
    ind = find(adj.roiWhole);
    FA{i}(ind)=atan2( abs(MT{i}(ind)),MZ{i}(ind) )/pi*180;
    
    display_bloch_simulation(MT{i},MZ{i},FA{i},BETA2{i},opt,adj,'LFA');
    
    errorLTA{i} = calculateError(MT{i},FA{i},BETA2{i},adj,opt, i);
end

% power & SAR
[avgp, maxp] = display_optimization_summary(rf,adj,spokes,sar,system);

%% PRINT

errorLTA = errorLTA{1};
if(opt.avg_target == 1)
    if(opt.refoc_pulse == 1)
        fprintf('\n ____________ BC Avg Refocusing________________\n')
    else
        fprintf('\n ____________ BC Avg Excitation________________\n')
    end
else
    if(opt.refoc_pulse == 1)
        fprintf('\n ____________ BC Center Refocusing________________\n')
    else
        fprintf('\n ____________ BC Center Excitation________________\n')
    end
end
fprintf('\n RMMSE \n')
fprintf('FA: %2.2f %%\n', 100*errorLTA.FArmmse);
fprintf('MT: %2.2f %%\n', 100*errorLTA.MTrmmse);
fprintf('beta2: %2.2f %%\n', 100*errorLTA.betarmmse);
fprintf('\n RMSE \n')
fprintf('MT: %2.2f %%\n', 100*errorLTA.MTrmse);
fprintf('beta2: %2.2f %%\n', 100*errorLTA.betarmse);

%% SAVE RESULTS

pulse.LFA.xrf=xrf;
pulse.LFA.avgp = avgp;
pulse.LFA.maxp = maxp;
pulse.LFA.rf_full=rf_full;
pulse.LFA.MT=MT;
pulse.LFA.MZ=MZ;
pulse.LFA.FA=FA;
pulse.LFA.BETA2=BETA2;
pulse.LFA.comptime=0;
systemm = system;
save('all', 'pulse', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA','sar');


end

