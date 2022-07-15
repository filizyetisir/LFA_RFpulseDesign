

function gpuStruct=prepare_gpu(adj,system,spokes)

gpuStruct=[];

if parallel.gpu.GPUDevice.isAvailable
    
%     % grab GPU
%     gpuStruct.gpuDeviceCount = gpuDeviceCount;
%     if gpuStruct.gpuDeviceCount > 1
%         fprintf('Found %d GPUs, please select the GPU to use:\n',gpuStruct.gpuDeviceCount)
%         for i=1:gpuStruct.gpuDeviceCount
%             g=gpuDevice(i);
%             fprintf('\t[%d]: ''%s'' has compute capability %s, %d threads per block, %d SMs\n', ...
%                 i,g.Name,g.ComputeCapability,g.MaxThreadsPerBlock,g.MultiprocessorCount);
%         end
%         n=input(sprintf('Enter GPU number: '));
%         gpuStruct.g=gpuDevice(n);
%     end

    gpuStruct.g=gpuDevice;
    
    % prepare kernel -- objective function and gradient
    fprintf('Preparing objective function kernel ...\n');
    gpuStruct.f0_LFA_gpuKernel = parallel.gpu.CUDAKernel('f0_LFA_cudaKernel_v3.ptx','f0_LFA_cudaKernel_v3.cu');
    
    M=128;
    N=ceil( adj.nnonzeropixels / M );
    gpuStruct.f0_LFA_gpuKernel.ThreadBlockSize = M;
    gpuStruct.f0_LFA_gpuKernel.GridSize = N;
        
    gpuStruct.a_re=gpuArray( ones(adj.nnonzeropixels,1) );
    gpuStruct.a_im=gpuArray( zeros(adj.nnonzeropixels,1) );
    gpuStruct.b_re=gpuArray( zeros(adj.nnonzeropixels,1) );
    gpuStruct.b_im=gpuArray( zeros(adj.nnonzeropixels,1) );
    
    gpuStruct.datot_dre_re = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.datot_dre_im = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.datot_dim_re = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.datot_dim_im = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    
    gpuStruct.dbtot_dre_re = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.dbtot_dre_im = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.dbtot_dim_re = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.dbtot_dim_im = gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    
    gpuStruct.a_forw_re=gpuArray( zeros(adj.nnonzeropixels,spokes.ntimes) );
    gpuStruct.a_forw_im=gpuArray( zeros(adj.nnonzeropixels,spokes.ntimes) );
    
    gpuStruct.b_forw_re=gpuArray( zeros(adj.nnonzeropixels,spokes.ntimes) );
    gpuStruct.b_forw_im=gpuArray( zeros(adj.nnonzeropixels,spokes.ntimes) );    

    % prepare kernel -- hessian
    fprintf('Preparing Hessian kernel ...\n');
    gpuStruct.hessian_LFA_noB0_cudaKernel = parallel.gpu.CUDAKernel('hessian_LFA_noB0_cudaKernel.ptx','hessian_LFA_noB0_cudaKernel.cu');
    
    M=128;
    N=ceil( adj.nnonzeropixels / M );
    gpuStruct.hessian_LFA_noB0_cudaKernel.ThreadBlockSize = M;
    gpuStruct.hessian_LFA_noB0_cudaKernel.GridSize = N;
    
    ncunk=system.ncoils * spokes.nspokes;
    
    gpuStruct.ha_re=gpuArray( zeros(adj.nnonzeropixels,2*ncunk,2*ncunk) );
    gpuStruct.ha_im=gpuArray( zeros(adj.nnonzeropixels,2*ncunk,2*ncunk) );
    
    gpuStruct.hb_re=gpuArray( zeros(adj.nnonzeropixels,2*ncunk,2*ncunk) );
    gpuStruct.hb_im=gpuArray( zeros(adj.nnonzeropixels,2*ncunk,2*ncunk) );

    gpuStruct.ja_re_re=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.ja_re_im=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.ja_im_re=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.ja_im_im=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    
    gpuStruct.jb_re_re=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.jb_re_im=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.jb_im_re=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    gpuStruct.jb_im_im=gpuArray( zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes) );
    
end















