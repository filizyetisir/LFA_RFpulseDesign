


function spokes=opt_gradients(spokes,system,opt,sar,adj,do_mls,SP_mean)
% function spokes=opt_gradients(spokes,system,opt,sar,adj,do_mls,SP_mean)

global A;
global spcoords;

spcoords=inf( 2*spokes.nspokes,1 );  % initial the sp coord vector with garbage to force computation of the system matrix
ncunk=system.ncoils*spokes.nspokes;

% START FROM LS/MLS PULSE
spokes=prepare_pulse(spokes,opt,system);
A=compute_system_matrix(adj,system,spokes,0,SP_mean);

optim_options=optimset('Display','none','MaxFunEvals',10^5,'TolFun',1e-4,'TolCon',1e-4,'MaxIter',5000, ...
    'Algorithm','interior-point','GradObj','on','GradConstr','on');

optim_options.Hessian='user-supplied';
optim_options.HessFcn=@(xrf,lambda)dualhess_STA(xrf,lambda,A,adj,spokes,system,sar,1,do_mls);

xrf0=zeros(2*ncunk,1);
[xrf fval exitflag output lambda]=fmincon(@(xrf)f0_STA(xrf,A,adj,do_mls),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,adj,spokes,sar,system),optim_options);

xsp0=zeros(2*spokes.nspokes,1);
for i=1:spokes.nspokes
    xsp0( (i-1)*2+1,1 )=spokes.spcoords(i,1);
    xsp0( (i-1)*2+2,1 )=spokes.spcoords(i,2);
end

% fac=max(abs(xrf));
fac=100;  % this makes the spoke locations on the same order than the RF values
x0=[xsp0*fac;xrf];

% JOINTLY OPTIMZIZE SPOKE POSITIONS AND RF

% test_grad(@(x)f0_STA_SP(x,adj,opt,spokes,system,do_mls,fac,SP_mean),rand(size(x0))*100)
% lambda=zeros(10000,1);
% test_hessian(@(x)f0_STA_SP(x,adj,opt,spokes,system,do_mls,fac,SP_mean),@(x)dualhess_STA_SP(x,lambda,fac,adj,spokes,sar,system,opt,SP_mean,do_mls,2),rand(size(x0))*100);

optim_options=optimset('Display','iter-detailed','MaxFunEvals',10^5,'TolFun',1e-4,'TolCon',1e-4,'TolX',1e-10, ...
    'MaxIter',20,'Algorithm','interior-point','GradObj','on','GradConstr','on');

% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(x,lambda)dualhess_STA_SP(x,lambda,fac,adj,spokes,sar,system,opt,SP_mean,do_mls,1);
options.Hessian='';

[x fval exitflag output lambda]=fmincon( @(x)f0_STA_SP(x,adj,opt,spokes,system,do_mls,fac,SP_mean),x0,[],[],[],[],[],[],@(x)fn_SP(x,adj,spokes,sar,system),optim_options );

for i=1:spokes.nspokes
    spokes.spcoords(i,1)=x( (i-1)*2+1,1 )/fac;
    spokes.spcoords(i,2)=x( (i-1)*2+2,1 )/fac;
end










