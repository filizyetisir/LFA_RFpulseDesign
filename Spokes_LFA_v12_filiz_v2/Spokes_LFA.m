function Spokes_LFA(filepath,opt,lims,i)
% function Spokes_LFA(filepath,opt)

global gpuStruct;


% Read parameters
fprintf('Reading pulse parameters ...\n');
[adj, system, opt, spokes, sar]=read_spokes_def_file(filepath,opt);
spokes.SS_grad_down = opt.gDown/100;

sar.lsarmax = lims.lsarlim;
sar.gsarmax = lims.gsarlim;
sar.ppavmax = lims.avpowlim;
system.umax = lims.rf_max;

% Make b0 zero for test purposes
if(opt.b0_zero == 1)
    adj.b0mapWhole = adj.b0mapWhole*0;
    adj.b0mapInner = adj.b0mapInner*0;
    adj.b0mapSpecial = adj.b0mapSpecial*0;
end



opt.nfreqs=size(opt.freqs,2);
spokes.verse_factor=opt.verse_factor;


% Prepare optimization
fprintf('Preparing optimization ...\n');
adj=prepare_optimization(opt,spokes,system,adj);
spokes=prepare_pulse(spokes,opt,system);  % prepare pulse on the system raster time (used to export the ini files)
spokes.q_gblips=get_Q_gradient_blips(spokes,adj,system);  % get rotation spinors (diagonal) for gradient blips


fprintf('Subpulse duration: %1.2f, total pulse duration: %1.2f\n', spokes.ntimesperspoke, spokes.ntimes)


% Read sar matrices
fprintf('Reading SAR matrices ...\n');
sar=prepare_sar(sar,spokes,system);


% Prepare GPU
if opt.use_gpu~=0
    gpuStruct=prepare_gpu(adj,system,spokes);
end



% Optimization options
tolfun = 1e-3;
tolx = 1e-6;
tolc = 1e-6;
optim_options=optimset('display','iter-detailed','MaxFunEvals',10^5,'TolFun',tolfun,'TolCon',tolc, 'TolX', tolx,'MaxIter',500,...
    'Algorithm','interior-point','GradObj','on','GradConstr','on');



if(spokes.nspokes == 1 && opt.use_bc_init_LFA == 1)
    [ind1, ind2] = size(adj.roiWhole);
    fil1 = opt.fil1; fil2 = opt.fil2;
    % give a 2 cm range
    range = round(0.01/adj.dx);
    ind1 = round(ind1/2)+fil1; ind2 = round(ind2/2)+fil2;
    
    bcmode = zeros(size(squeeze(adj.b1maps(1,:,:))));
    %pha = [0.78; 1.57 ;2.35 ;3.14; 3.93; 4.71 ;5.50; 6.28];
    for c = 1:size(adj.b1maps,1)
        bcmode = bcmode + squeeze(adj.b1maps(c,:,:).*exp(-1j*angle(adj.b1maps(c,ind1,ind2))));
        %bcmode = bcmode + squeeze(adj.b1maps(c,:,:));%*exp(-1j*(pha(c)));
    end
    
    % Find the scaling factor to achieve the target
    gamma = 42.576e6;
    ind1_range = [ind1-range:ind1+range];
    ind2_range = [ind2-range:ind2+range];
    tmp = abs(bcmode(ind1_range,ind2_range));
    tmp = mean(tmp(:));
    scaleFac = opt.tfa/(360*system.deltat*gamma*sum(spokes.sinc)*tmp);
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
    xrf_STA = [real(rf); imag(rf)];
    
else

    
% STA design
fprintf('\n*********  S T A   P U L S E   D E S I G N  *********\n');
[xrf_STA, opt, ncunk, rf_STA_full, MT_STA, MZ_STA, FA_STA, BETA2_STA, comptime_STA, spokes, errorSTA, avgp_STA, maxp_STA]=design_STA(system,opt,adj,spokes,sar,optim_options);

end

% Calculate SAR
sar_STA = calcVOPSAR(sar,spokes,system,xrf_STA);




fprintf('\n*********  L F A   R E F I N E M E N E T  ( W I T H   B 0 )  *********');
[xrf_LFA, rf_LFA_full, MT_LFA, MZ_LFA, FA_LFA, BETA2_LFA, timeAll, optimality, errorLTA, MTtimely, avgp_LTA, maxp_LTA]=design_LFA(xrf_STA,adj,system,spokes,opt,sar,optim_options);

% Calculate SAR
sar_LTA = calcVOPSAR(sar,spokes,system,xrf_LFA);


sar.STA = sar_STA;
sar.LTA = sar_LTA;

% Save Results
pulse.g=spokes.g;
pulse.k=spokes.k;

if(spokes.nspokes == 1 && opt.use_bc_init_LFA == 1)
    pulse.STA.xrf=xrf_LFA;
    pulse.STA.rf_full=rf_LFA_full;
    pulse.STA.MT=MT_LFA;
    pulse.STA.MZ=MZ_LFA;
    pulse.STA.FA=FA_LFA;
    pulse.STA.BETA2=BETA2_LFA;
    pulse.STA.comptime=0;
    errorSTA = errorLTA;
else
    pulse.STA.xrf=xrf_STA;
    pulse.STA.rf_full=rf_STA_full;
    pulse.STA.MT=MT_STA;
    pulse.STA.MZ=MZ_STA;
    pulse.STA.FA=FA_STA;
    pulse.STA.BETA2=BETA2_STA;
    pulse.STA.comptime=comptime_STA;
end

pulse.LFA.xrf=xrf_LFA;
pulse.LFA.rf_full=rf_LFA_full;
pulse.LFA.MT=MT_LFA;
pulse.LFA.MZ=MZ_LFA;
pulse.LFA.FA=FA_LFA;
pulse.LFA.BETA2=BETA2_LFA;
pulse.LFA.timeAll=timeAll;
pulse.LFA.optimality = optimality;
pulse.LFA.avgp = avgp_LTA;
pulse.LFA.maxp = maxp_LTA;
systemm = system;


save('all', 'pulse', 'sar', 'opt', 'adj', 'systemm', 'spokes', 'errorLTA', 'errorSTA', 'MTtimely');













