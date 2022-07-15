

function [xrf rf_full MT MZ FA BETA2 comptime]=design_LFA_noB0(xrf0,adj,system_opt,system,spokes_opt,spokes,opt,sar,optim_options,refoc_pulse)
% function [xrf rf_full MT MZ FA BETA2 comptime]=design_LFA_noB0(xrf0,adj,system_opt,system,spokes_opt,spokes,opt,sar,optim_options,refoc_pulse)

% lambda=zeros(1000,1);
% test_grad( @(xrf)f0_LFA_noB0(xrf,adj,spokes_opt,system_opt,opt),xrf0 );
% test_hessian( @(xrf)f0_LFA_noB0(xrf,adj,spokes_opt,system_opt,opt),@(xrf)dualhess_LFA_noB0(xrf,lambda,adj,system_opt,spokes_opt,opt,sar,2),xrf0 );
% test_hessian( @(xrf)f0_LFA_noB0(xrf,adj,spokes_opt,system_opt,opt),@(xrf)hessian_LFA_noB0_mexInterf(xrf,lambda,adj,system_opt,spokes_opt,opt,sar,2),xrf0 );

optim_options.Hessian='user-supplied';
optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexInterf(xrf,lambda,adj,system_opt,spokes_opt,opt,sar,1);
% optim_options.Hessian='';

tic;
[xrf fval exitflag output lambda]=fmincon(@(xrf)f0_LFA_noB0(xrf,adj,spokes_opt,system_opt,opt),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,adj,spokes_opt,sar,system_opt),optim_options);
comptime=toc;

% smal adjustment due to posisble small difference between the G plateau of
% both spokes structures (should be pretty close though...)
xrf=xrf*(spokes.Gss/spokes_opt.Gss);

% save pulse (using the system raster time, not the optimization raster
% time)
ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk,1) + 1j*xrf(ncunk+1:2*ncunk,1);

% rf_full=write_and_display_pulse(rf,opt,spokes_opt,system_opt,'LFA -- no B0');
rf_full=write_and_display_pulse(rf,opt,spokes,system,'LFA -- no B0');

if exist('pTXArbitrary.ini','file')
    movefile('pTXArbitrary.ini','pTXArbitrary_LTA_noB0.ini');
else
    warning('The LFA (no B0) pulse file pTXArbitrary.ini does not exist.');
end


% Bloch sim
if opt.refoc_pulse==1
    M0=[0;0;1];
else
    M0=[0;0;1];
end

% [mx my mz BETA2]=bloch_simulation(rf_full,adj,spokes_opt,system_opt,0,M0);
[mx my mz BETA2]=bloch_simulation(rf_full,adj,spokes,system,0,M0);

MT=mx + 1j*my;
MZ=mz;
FA=zeros(size(adj.roi));
% FA(adj.ind)=acos(mz(adj.ind))/pi*180;    
FA(adj.ind)=atan2( abs(MT(adj.ind)),MZ(adj.ind) )/pi*180;
display_bloch_simulation(MT,MZ,FA,BETA2,opt,'LFA -- no B0');

% power & SAR
% display_optimization_summary(rf,adj,spokes_opt,sar,system_opt);
display_optimization_summary(rf,adj,spokes,sar,system);















