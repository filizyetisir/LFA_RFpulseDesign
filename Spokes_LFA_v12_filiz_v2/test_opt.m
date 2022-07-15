%test for the optimization methods
clear;
load data;%A opt sar spokes system

options = optimset('Largescale','off','Display','iter-detailed','MaxFunEvals',10^5,'TolFun',1e-6,'TolCon',1e-6,'MaxIter',5000,...
    'Algorithm','interior-point','Hessian','user-supplied','HessFcn',@(xrf,lambda)dualhess(xrf,lambda,A,opt,spokes,sar),...
    'GradObj','on','GradConstr','on','SubproblemAlgorithm','ldl-factorization');                                    %what's this...
ncunk=opt.nrows*opt.ncoilspr*spokes.nspokes;
xrf0=zeros(2*ncunk,1);
[xrf fval exitflag output lambda]=fmincon(@(xrf)f0(xrf,A,opt),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,opt,spokes,sar),options);

% process results
rf=xrf(1:ncunk,1) + 1j*xrf(ncunk+1:2*ncunk,1);
display_optimization_summary(rf,opt,spokes,sar);
rf2=write_and_display_pulse(rf,opt,spokes,system);
write_and_display_magnetization(rf,opt,system,spokes);
[mx my mz]=bloch_simulation(rf2,opt,spokes,system);
write_and_display_bloch_simulation(mx,my,mz,opt);
