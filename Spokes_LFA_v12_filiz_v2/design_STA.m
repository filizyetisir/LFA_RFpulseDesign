

function [xrf, opt, ncunk, rf_full, MT, MZ, FA, BETA2, comptime, spokes, error, avgp, maxp]=design_STA(system,opt,adj,spokes,sar,optim_options)
% function [xrf opt ncunk rf_full MT MZ FA BETA2 comptime spokes]=design_STA(system,opt,adj,spokes,sar,optim_options)

% scale SAR limit to ensure feasibility even after scaling
% of the STA pulse
if opt.refoc_pulse
    tfa=180.0;
else
    tfa=opt.tfa;
end   
FA_FACTOR=30/tfa;  % ratio of initial STA FA and target LFA FA
add_SF = 1;
system.umax=system.umax * FA_FACTOR*add_SF;
sar.gsarmax=sar.gsarmax * FA_FACTOR^2*(add_SF^2);
sar.lsarmax=sar.lsarmax * FA_FACTOR^2*(add_SF^2);
sar.ppavmax = sar.ppavmax * FA_FACTOR^2*(add_SF^2);

% OPTIMIZATION SPOKE POSITIONS IF NEEDED

if opt.opt_grad==1    
    fprintf('\n*********  O P T I M I Z I N G   S P O K E S   P O S I T I O N S  *********\n');
    display=opt.display;
    opt.display='off';    
    
    spokes=opt_gradients_gridsearch(spokes,system,opt,sar,adj,opt.do_mls,1,10);    
    
    spokes=opt_gradients(spokes,system,opt,sar,adj,opt.do_mls,1);    
    
    spokes=prepare_pulse(spokes,opt,system);
    
    opt.display=display;
end

ncunk=spokes.nspokes * system.ncoils;
opt.A_STA=zeros(opt.nfreqs*adj.nnonzeropixels,ncunk);

%adj.target_STA=repmat(adj.target_STA,[opt.nfreqs 1]);

tic
for i=1:opt.nfreqs
    ind=(i-1)*adj.nnonzeropixels + 1 : i*adj.nnonzeropixels;
    opt.A_STA(ind,:)=compute_system_matrix(adj,system,spokes,opt.freqs(i),1) * sqrt(opt.freq_weights(i));
    
    adj.target_STA(ind)=adj.target_STA(ind) * sqrt(opt.freq_weights(i));
end
comptime1 = toc;
fprintf('A matrix computation time: %2.2f\n',comptime1)

xrf0=zeros(2*ncunk,1);

if(spokes.nspokes == 1 && opt.use_bc_init_STA == 1)
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
    scaleFac = 30/(360*system.deltat*gamma*sum(spokes.sinc)*tmp);
    if(scaleFac > system.umax)
        scaleFac = system.umax;
        %error('BC mode RF pulse exceeds maximum voltage limit. Make the pulse duration longer')
    end
    
    % Form the BC mode pulse
    pha = angle(adj.b1maps(:,ceil(end/2+fil1),ceil(end/2+fil2)));
    if(opt.refoc_pulse == 0)
        rf = scaleFac*ones(8,1).*exp(-1j*pha);
    else
        rf = scaleFac*ones(8,1).*exp(-1j*pi/2).*exp(-1j*(pha));
    end
    xrf0 = [real(rf); imag(rf)];
end

optim_options.Hessian='user-supplied';
optim_options.HessFcn=@(xrf,lambda)dualhess_STA(xrf,lambda,opt.A_STA,adj,spokes,system,sar,1,opt.do_mls,opt);    
%optim_options.HessFcn=@(xrf,lambda)dualhess_STA_withAvgPower(xrf,lambda,opt.A_STA,adj,spokes,system,sar,1,opt.do_mls,opt);    
%optim_options.Hessian='';
optim_options.GradObj = 'on';


%% TESTING

% lambda=zeros(1000,1);
% xrf0 = rand(32,1)*150; % for 2 spokes

% test_grad( @(x)f0_STA(x,opt.A_STA,adj,opt.do_mls),xrf0 );

% test_hessian only tets for the hessian calculation of the oibjective
% function hence we make lambda(lagrange multipliers) all zeros, 
% which is related to the constraints 
% test_hessian( @(x)f0_STA(x,opt.A_STA,adj,opt.do_mls),@(x)dualhess_STA(x,lambda,opt.A_STA,adj,spokes,system,sar,0,opt.do_mls,opt),xrf0 );

% test_gradient_fn( @(x)fn(x,adj,spokes,sar,system,opt), xrf0)

%% 

tic
[xrf, fval, exitflag, output, lambda]=fmincon(@(xrf)f0_STA(xrf,opt.A_STA,adj,opt.do_mls),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,adj,spokes,sar,system,opt),optim_options);
comptime=toc;
fprintf('\nSTA design time: %2.2f\n',comptime)

% scale STA design to LFA
xrf=xrf/FA_FACTOR;

% smal adjustment due to posisble small difference between the G plateau of
% both spokes structures (should be pretty close though...)
xrf=xrf*(spokes.Gss/spokes.Gss);

% save pulse
rf=xrf(1:ncunk,1) + 1j*xrf(ncunk+1:2*ncunk,1);
rf_full=write_and_display_pulse(rf,opt,spokes,system,'scaled STA');

if exist('pTXArbitrary.ini','file')
    movefile('pTXArbitrary.ini','pTXArbitrary_STA.ini');
else
    warning('The STA pulse file pTXArbitrary.ini does not exist.');
end




nsims=length(opt.freqs);  % must be odd to simulate b0_off=0 Hz (Larmor)
% b0_offsets=linspace(-100,100,nsims)';
b0_offsets=opt.freqs;

MT=cell(nsims,1);
MTtimely=cell(nsims,1);
MZ=cell(nsims,1);
FA=cell(nsims,1);
BETA2=cell(nsims,1);
rmses=zeros(nsims,1);

% Bloch sim
if opt.refoc_pulse
    M0=[0;0;1];
else
    M0=[0;0;1];
end
% [mx my mz BETA2]=bloch_simulation(rf_full,adj,spokes,system,0,M0);
% 
% 
% MT=mx + 1j*my;
% MZ=mz;
% FA=zeros(size(adj.roi));
% ind = find(adj.roiWhole);   
% FA(ind)=atan2( abs(MT(ind)),MZ(ind) )/pi*180;
% display_bloch_simulation(MT,MZ,FA,BETA2,opt,adj,'STA');
% 
% error = calculateError(MT, FA, BETA2, adj, opt);


for i=1:size(b0_offsets,2)

    [mx, my, mz, BETA2{i}, mxytimely]=bloch_simulation(rf_full,adj,spokes,system,b0_offsets(i),M0);

    MT{i}=mx + 1j*my;
    MZ{i}=mz;
    
    MTtimely{i} = mxytimely;
    
    FA{i}=zeros(size(adj.roi));  
    ind = find(adj.roiWhole);   
    FA{i}(ind)=atan2( abs(MT{i}(ind)),MZ{i}(ind) )/pi*180;
       
    display_bloch_simulation(MT{i},MZ{i},FA{i},BETA2{i},opt,adj,'STA');
    
    error{i} = calculateError(MT{i}, FA{i}, BETA2{i}, adj, opt, i);
end





% power & SAR
[avgp, maxp] = display_optimization_summary(rf,adj,spokes,sar,system);


