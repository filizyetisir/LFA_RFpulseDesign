

function [xrf, rf_full, MT, MZ, FA, BETA2, timeAll, optimality, error, MTtimely, avgp, maxp]=design_LFA(xrf0,adj,system,spokes,opt,sar,optim_options)
%function [xrf, rf_full, MT, MZ, FA, BETA2, comptime, error, MTtimely]=design_LFA(xrf0,adj,system,spokes,opt,sar,optim_options)

% lambda=ones(1000,1);
% test_grad( @(xrf)f0_LFA(xrf,adj,spokes,system,opt),xrf0 );
% test_grad( @(xrf)f0_LFA_mexCudaInterf(xrf,adj,spokes,system,opt),xrf0 );
% test_hessian( @(xrf)f0_LFA_mexCudaInterf(xrf,adj,spokes,system,opt),@(xrf)dualhess_LFA_noB0(xrf,lambda,adj,system,spokes,opt,sar,0),xrf0 );
% test_hessian( @(xrf)f0_LFA_mexCudaInterf(xrf,adj,spokes,system,opt),@(xrf)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,0),xrf0 );

if(opt.useH == 1)
    optim_options.Hessian='user-supplied';
    optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);  
else
    optim_options.Hessian='';
end
if(opt.useJ == 1)
    optim_options.GradObj = 'on';
else
    optim_options.GradObj = 'off';
end

if(opt.refoc_pulse == 1)
    tolf = 1e-2;
else
    tolf = 1e-1;
end
optim_options.MaxIter=500;
optim_options.TolFun=tolf;
optim_options.TolCon=1e-6;
optim_options.TolX = 1e-6;

% This is for pulse design time and optimality plot for the paper
global firstorderopt;
firstorderopt = [];
global time;
time = [];
optim_options.OutputFcn = @outfun;


%% TESTING

% lambda=zeros(1000,1);
% xrf0 = rand(32,1)*150; % for 2 spokes

% test_grad( @(x)f0_LFA_mexCudaInterf(x,adj,spokes,system,opt),xrf0 );

% test_hessian only tets for the hessian calculation of the oibjective
% function hence we make lambda(lagrange multipliers) all zeros, 
% which is related to the constraints 
% test_hessian( @(x)f0_LFA_mexCudaInterf(x,adj,spokes,system,opt),@(x)hessian_LFA_noB0_mexCudaInterf(x,lambda,adj,system,spokes,opt,sar,0),xrf0 );

% test_gradient_fn( @(x)fn(x,adj,spokes,sar,system,opt), xrf0)

%% 

tic
[xrf, fval, exitflag, output, lambda]=fmincon(@(xrf)f0_LFA_mexCudaInterf(xrf,adj,spokes,system,opt),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,adj,spokes,sar,system,opt),optim_options);
comptimeLTA=toc;
fprintf('LTA design time1: %2.2f\n',comptimeLTA)

% This is for pulse design time and optimality plot for the paper
timeAll(1) = time(1); 
for i = 2:length(time)
    timeAll(i) = timeAll(i-1) + time(i); 
end
comptime = timeAll;
fprintf('LTA design time2: %2.2f\n',comptime(end))
optimality = firstorderopt;


% small adjustment due to posisble small difference between the G plateau of
% both spokes structures (should be pretty close though...)
xrf=xrf*(spokes.Gss/spokes.Gss);

% save pulse
ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk,1) + 1j*xrf(ncunk+1:2*ncunk,1);
rf_full=write_and_display_pulse(rf,opt,spokes,system,'LFA');

if exist('pTXArbitrary.ini','file')
    movefile('pTXArbitrary.ini','pTXArbitrary_LTA.ini');
else
    warning('The LFA pulse file pTXArbitrary.ini does not exist.');
end

% Bloch sim
if opt.refoc_pulse==1
    M0=[0;0;1];
else
    M0=[0;0;1];
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
for i=1:size(b0_offsets,2)

    [mx, my, mz, BETA2{i}, mxytimely]=bloch_simulation(rf_full,adj,spokes,system,b0_offsets(i),M0);

    MT{i}=mx + 1j*my;
    MZ{i}=mz;
    
    MTtimely{i} = mxytimely;
    
    FA{i}=zeros(size(adj.roi));  
    ind = find(adj.roiWhole);   
    FA{i}(ind)=atan2( abs(MT{i}(ind)),MZ{i}(ind) )/pi*180;
       
    display_bloch_simulation(MT{i},MZ{i},FA{i},BETA2{i},opt,adj,'LFA');
    
    error{i} = calculateError(MT{i}, FA{i}, BETA2{i}, adj, opt, i);
end


% power & SAR
[avgp, maxp] = display_optimization_summary(rf,adj,spokes,sar,system);


















