
% You can make the tolerance for objective function, constraints and the
% step size of x 10^-2. It has no difference from when its 10^-4 in terms
% of the rmmse and rmse resulting for both excitation and refocusing pulses
% designed with LS and MLS approaches. 

% For Jacobian, Hessian, GPU
% in runSpokes
% usegpu = 1 
% in f0_LFA_mexCudaInterf
% make sure there is both f and df in the output and comp_grad = 1 inside 
% In design_LFA
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);   
% optim_options.GradObj = 'on';

% For Jacobian, Hessian, no GPU
% clear all
% in runSpokes
% use_gpu = 0
% in f0_LFA_mexCudaInterf
% make sure there is both f and df in the output and comp_grad = 1 inside 
% in design_LFA
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);   
% optim_options.GradObj = 'on';

% For Jacobian, no Hessian, GPU
% Same as above except
% in runSpokes
% use_gpu = 1;
% in f0_LFA_mexCudaInterf
% make sure there is both f and df in the output and comp_grad = 1 inside 
% in design_LFA
% optim_options.GradObj = 'on';
% Comment out 
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);  
% Make 
% optim_options.Hessian = '';

% For Jacobian, no Hessian, no GPU
% clear all
% in run_Spokes
% use_gpu = 0;
% in f0_LFA_mexCudaInterf
% make sure there is both f and df in the output and comp_grad = 1 inside 
% in design_LFA
% optim_options.GradObj = 'on';
% Comment out 
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);  
% Make 
% optim_options.Hessian = '';

% For no Jacobian, no Hessian, no GPU
% clear all
% in run_Spokes
% use_gpu = 0;
% in f0_LFA_mexCuda_Interf
% delete df, make comp_grad = 0;
% in design_LFA
% make optim_options.GradObj = 'off';
% Comment out 
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);  
% Make 
% optim_options.Hessian = '';

% For no Jacobian, no Hessian, GPU
% clear all
% in run_Spokes
% use_gpu = 1;
% in f0_LFA_mexCuda_Interf
% delete df, make comp_grad = 0;
% in design_LFA
% make optim_options.GradObj = 'off';
% Comment out 
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);  
% Make 
% optim_options.Hessian = '';


% BELOW two are not possible, you need jacbian for Hessian

% For no Jacobian, Hessian, no GPU
% clear all
% in run_Spokes
% use_gpu = 0;
% in f0_LFA_mexCuda_Interf
% delete df, make comp_grad = 0;
% in design_LFA
% make optim_options.GradObj = 'off';
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);   

% For no Jacobian, Hessian, GPU
% clear all
% in run_Spokes
% use_gpu = 1;
% in f0_LFA_mexCuda_Interf
% delete df, make comp_grad = 0;
% in design_LFA
% make optim_options.GradObj = 'off';
% optim_options.Hessian='user-supplied';
% optim_options.HessFcn=@(xrf,lambda)hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,1);   



