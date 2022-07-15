

% read matrices from dump file
mfilepath='Y:\HFSS\StackZLoopCoils\1row_8coils_2\PPD_decoupled_trans100\performance_pd_schur_1290_vops\matrices_dump.txt'
[npixels nspokes nchannels nconst x0 lambdas0 mtarget sm const f0hess]=read_matrices_dump(mfilepath);

% optimize
tic
options = optimset('Largescale','off','Display','iter','MaxFunEvals',10^5,'TolFun',10^-10,'TolCon',10^-10,'MaxIter',5000,...
    'Algorithm','interior-point','Hessian','user-supplied','HessFcn',@(x,lambda)dualhess(x,lambda,nchannels,nspokes,f0hess,nconst,const),...
    'GradObj','on','GradConstr','on','SubproblemAlgorithm','ldl-factorization');
[x fval exitflag output lambda]=fmincon(@(x)f0(sm,mtarget,x),x0,[],[],[],[],[],[],@(x)fn(x,nspokes,nconst,nchannels,const,10,1),options);
toc
% x
% lambda.ineqnonlin

