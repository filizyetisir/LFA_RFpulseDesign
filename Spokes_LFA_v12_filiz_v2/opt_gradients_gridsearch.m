


function spokes=opt_gradients_gridsearch(spokes,system,opt,sar,adj,do_mls,SP_mean,npoints)
% function spokes=opt_gradients_gridsearch(spokes,system,opt,sar,adj,do_mls,SP_mean,npoints)


spokes_=spokes;

lx=max(adj.nonzeropixels(:,4)) - min(adj.nonzeropixels(:,4));
ly=max(adj.nonzeropixels(:,5)) - min(adj.nonzeropixels(:,5));

grid_kx=linspace(-1.0/lx,1.0/lx,npoints);
grid_ky=linspace(-1.0/ly,1.0/ly,npoints);

optim_options=optimset('Display','none','MaxFunEvals',10^5,'TolFun',1e-4,'TolCon',1e-4,'MaxIter',5000, ...
    'Algorithm','interior-point','GradObj','on','GradConstr','on');


fprintf('SPOKE LOCATION GRID SEARCH ');
sbuf=[];
for sp=2:spokes.nspokes  % position one spoke at a time
    
    spokes_.nspokes=sp;
    spokes_.spcoords = spokes.spcoords(end-sp+1:end,:);
    
    xrf0=zeros(2*system.ncoils*sp,1);

    % loop over grid points
    f0_grid=zeros(npoints);
    f0=inf;
    fval=inf;
    for i=1:npoints
        for j=1:npoints    
                
            spokes_.spcoords(end-sp+1,:)=[ grid_kx(i) grid_ky(j) 0.0 ];   
            
            spokes_=prepare_pulse(spokes_,opt,system);
            A=compute_system_matrix(adj,system,spokes_,0,SP_mean);
            
            optim_options.Hessian='user-supplied';
            optim_options.HessFcn=@(xrf,lambda)dualhess_STA(xrf,lambda,A,adj,spokes_,system,sar,1,do_mls);
            
            [xrf fval exitflag output lambda]=fmincon(@(xrf)f0_STA(xrf,A,adj,do_mls),xrf0,[],[],[],[],[],[],@(xrf)fn(xrf,adj,spokes_,sar,system),optim_options);
            
            f0_grid(i,j)=fval;
            
            if fval<f0
                f0=fval;
                i_min=i;
                j_min=j;
            end
            
            for k=1:size(sbuf,2)
                fprintf('\b');
            end
            sbuf=sprintf('[%d/%d  f0=%f  fmin=%f]',(sp-2)*npoints^2 + (i-1)*npoints + j,(spokes.nspokes-1)*npoints^2,fval,f0);
            fprintf(sbuf);
            
        end
    end
    
    f=figure;
    surf(grid_kx,grid_ky,f0_grid); box off;
    xlabel('kx [1/m]'); ylabel('ky [1/m]'); title(sprintf('Grid search for spoke %d/%d',sp,spokes.nspokes));
    
    spokes.spcoords(end-sp+1,:)=[ grid_kx(i_min) grid_ky(j_min) 0.0 ];  
end
fprintf('\n');







