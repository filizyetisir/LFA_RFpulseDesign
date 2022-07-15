


function h=hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,solver)
% function h=hessian_LFA_noB0_mexCudaInterf(xrf,lambda,adj,system,spokes,opt,sar,solver)


global gpuStruct;

ncunk=size(xrf,1)/2;
nchannels=ncunk/spokes.nspokes;

rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);


% HESSIAN OF THE OBJECTIVE FUNCTION
btotspokes=zeros(adj.nnonzeropixels,spokes.nspokes);
for i=1:spokes.nspokes
    sp_ind=(i-1)*system.ncoils+1:i*system.ncoils;
    if isempty(gpuStruct)    
        btotspokes(:,i)=adj.b1maps_ * rf(sp_ind) * spokes.sumsinc; 
    else
        btotspokes(:,i)=gather( gpuArray(adj.b1maps_)*gpuArray(rf(sp_ind)) )*spokes.sumsinc;  % GPU!!
    end
end

if isempty(gpuStruct)  % no GPU --> use MEX implementation
    
    inds=1:adj.nnonzeropixels;
    
    opt_metric_tar = complex( real(adj.target_LFA(inds)),imag(adj.target_LFA(inds)) );
    h_OBJ=hessian_LFA_noB0_mex(btotspokes,spokes.q_gblips,adj.b1maps_,opt_metric_tar,spokes.sumsinc,system.deltat,opt.do_mls,opt.refoc_pulse);
    
else  % GPU
        
    % fast forward and backward Bloch sim
    [ gpuStruct.a_re,gpuStruct.a_im,gpuStruct.b_re,gpuStruct.b_im, ...
        gpuStruct.ha_re,gpuStruct.ha_im, ...
        gpuStruct.hb_re,gpuStruct.hb_im, ...
        gpuStruct.ja_re_re,gpuStruct.ja_re_im,gpuStruct.ja_im_re,gpuStruct.ja_im_im, ...
        gpuStruct.jb_re_re,gpuStruct.jb_re_im,gpuStruct.jb_im_re,gpuStruct.jb_im_im ] = feval( gpuStruct.hessian_LFA_noB0_cudaKernel, ...
            gpuStruct.a_re,gpuStruct.a_im,gpuStruct.b_re,gpuStruct.b_im, ...
            gpuStruct.ha_re,gpuStruct.ha_im, ...
            gpuStruct.hb_re,gpuStruct.hb_im, ...
            gpuStruct.ja_re_re,gpuStruct.ja_re_im,gpuStruct.ja_im_re,gpuStruct.ja_im_im, ...
            gpuStruct.jb_re_re,gpuStruct.jb_re_im,gpuStruct.jb_im_re,gpuStruct.jb_im_im, ...
            gpuStruct.a_forw_re,gpuStruct.a_forw_im,gpuStruct.b_forw_re,gpuStruct.b_forw_im, ...
            real(btotspokes),imag(btotspokes), ...
            real(spokes.q_gblips),imag(spokes.q_gblips), ...
            real(adj.b1maps_),imag(adj.b1maps_), ...
            adj.nnonzeropixels,spokes.nspokes,system.ncoils,spokes.sumsinc,system.deltat );
    
    a = gpuArray( gpuStruct.a_re + 1j*gpuStruct.a_im );
    b = gpuArray( gpuStruct.b_re + 1j*gpuStruct.b_im );
    
    ja_re = gpuArray( reshape( gpuStruct.ja_re_re + 1j*gpuStruct.ja_re_im,adj.nnonzeropixels,ncunk ) );
    ja_im = gpuArray( reshape( gpuStruct.ja_im_re + 1j*gpuStruct.ja_im_im,adj.nnonzeropixels,ncunk ) );
        
    jb_re = gpuArray( reshape( gpuStruct.jb_re_re + 1j*gpuStruct.jb_re_im,adj.nnonzeropixels,ncunk ) );
    jb_im = gpuArray( reshape( gpuStruct.jb_im_re + 1j*gpuStruct.jb_im_im,adj.nnonzeropixels,ncunk ) );

    ha = gpuArray( reshape( gpuStruct.ha_re + 1j*gpuStruct.ha_im,adj.nnonzeropixels,(2*ncunk)^2 ) );
    hb = gpuArray( reshape( gpuStruct.hb_re + 1j*gpuStruct.hb_im,adj.nnonzeropixels,(2*ncunk)^2 ) );

    if opt.refoc_pulse==0  % excitation pulse
        if opt.do_mls==1  % MLS
            opt_metric=abs(a).^2 - abs(b).^2;  % mz=|alpha|^2-|beta|^2
        else  % LS
            opt_metric=2.0*conj(a).*b;  % 2*conj(alpha)*beta
        end
    else  % refocusing pulse
        if opt.do_mls==1  % MLS
            opt_metric=abs(b).^2;  % |beta|^2
        else  % LS
            opt_metric=b.^2;  % beta^2
        end
    end
    
    inds=1:adj.nnonzeropixels;
    y=opt_metric - adj.target_LFA(inds);
    
    % Jacobian terms of total Hessian
    a_big=repmat(a,[1 ncunk]);
    b_big=repmat(b,[1 ncunk]);

    if opt.refoc_pulse==0  % excitation pulse
        if opt.do_mls==1  % MLS
            tmp1=2.0*real( ja_re.*conj(a_big) - jb_re.*conj(b_big) );
            tmp2=2.0*real( ja_im.*conj(a_big) - jb_im.*conj(b_big) );
        else  % LS
            tmp1=2.0*( conj(ja_re).*b_big + jb_re.*conj(a_big) );
            tmp2=2.0*( conj(ja_im).*b_big + jb_im.*conj(a_big) );
        end
    else  % refocusing pulse
        if opt.do_mls==1  % MLS
            tmp1=2.0*real( jb_re.*conj(b_big) );
            tmp2=2.0*real( jb_im.*conj(b_big) );
        else  % LS
            tmp1=2.0*jb_re.*b_big;
            tmp2=2.0*jb_im.*b_big;
        end
    end

    h11=gather( tmp1'*tmp1 );
    h12=gather( tmp1'*tmp2 );
    h22=gather( tmp2'*tmp2 );

    h1=2.0*real([ h11   h12;
                  h12.' h22 ] );
              
    % Hessian terms of total Hessian
    if opt.refoc_pulse==0  % excitation pulse
        if opt.do_mls==1  % MLS
            tmp1=reshape( (y.*a)'*ha - (y.*b)'*hb,2.0*ncunk,2.0*ncunk );
            tmp2=ja_re'*(repmat(y,[1 ncunk]).*ja_re) - jb_re'*(repmat(y,[1 ncunk]).*jb_re);
            tmp3=ja_im'*(repmat(y,[1 ncunk]).*ja_re) - jb_im'*(repmat(y,[1 ncunk]).*jb_re);
            tmp4=ja_im'*(repmat(y,[1 ncunk]).*ja_im) - jb_im'*(repmat(y,[1 ncunk]).*jb_im);
        else  % LS
            tmp1=reshape( (y.*conj(b))'*conj(ha) + (y.*a)'*hb,2.0*ncunk,2.0*ncunk );        
            tmp2=ja_re'*(repmat(conj(y),[1 ncunk]).*jb_re) + jb_re.'*(repmat(conj(y),[1 ncunk]).*conj(ja_re));
            tmp3=ja_im'*(repmat(conj(y),[1 ncunk]).*jb_re) + jb_im.'*(repmat(conj(y),[1 ncunk]).*conj(ja_re));
            tmp4=ja_im'*(repmat(conj(y),[1 ncunk]).*jb_im) + jb_im.'*(repmat(conj(y),[1 ncunk]).*conj(ja_im));        
        end
    else  % refocusing pulse
        if opt.do_mls==1  % MLS
            tmp1=reshape( (y.*b)'*hb,2.0*ncunk,2.0*ncunk );
            tmp2=jb_re'*(repmat(y,[1 ncunk]).*jb_re);
            tmp3=jb_im'*(repmat(y,[1 ncunk]).*jb_re);
            tmp4=jb_im'*(repmat(y,[1 ncunk]).*jb_im);
        else  % LS
            tmp1=reshape( (y.*conj(b))'*hb,2.0*ncunk,2.0*ncunk );
            tmp2=jb_re.'*(repmat(conj(y),[1 ncunk]).*jb_re);
            tmp3=jb_im.'*(repmat(conj(y),[1 ncunk]).*jb_re);
            tmp4=jb_im.'*(repmat(conj(y),[1 ncunk]).*jb_im);
        end
    end
    
    h2=4.0*real( gather( tmp1+[ tmp2 tmp3.';tmp3 tmp4 ] )  );

    % Hessian of objective function
    h_OBJ=h1+h2;
    
end

h_OBJ=h_OBJ*sum(opt.freq_weights,2);  % since we don't model B0 in the Hessian, the nfreqs designs have the same Hessian



% HESSIAN OF THE SAR AND POWER CONSTRAINTS
if solver==1  % Matlab solver
    LagMult=lambda.ineqnonlin;
else
    LagMult=lambda;
end

if(opt.sarcons)
    constant = spokes.sumsincsq/spokes.ntimes;
    % % Hessian of the constraints (SAR)
    subh=zeros(system.ncoils,system.ncoils);
    nconst=sar.nvop+1;
    subh=2*constant*squeeze( sum( sar.sarmats .* repmat(LagMult(1:nconst),[1 nchannels nchannels]),1 ) );
    ch=kron( eye(spokes.nspokes),subh ).*(system.ncoils);
    h_SAR=[real(ch) -imag(ch)';imag(ch)' real(ch)];
end

% Hessian of the constraints (max. power)
if(opt.sarcons)
    tmp=LagMult(nconst+1:nconst+ncunk);
else
    tmp=LagMult(1:ncunk);
end
tmp=repmat(tmp,[2 1]); % real and imaginary
h_MAX_POW=2.0*diag(tmp);


% Added 8/18/2020
% Hessian of the constraints (avg power)
if(opt.sarcons)
    tmp=LagMult(nconst+ncunk+1:nconst+ncunk+nchannels);
else
    tmp=LagMult(ncunk+1: ncunk+nchannels);
end
tmp=repmat(tmp,[2*spokes.nspokes 1]); % multiple spokes and real and imaginary
h_AV_POW=2.0*diag(tmp);


if(opt.sarcons)
    h=h_OBJ + h_SAR + h_MAX_POW + h_AV_POW;
else
    h=h_OBJ + h_MAX_POW + h_AV_POW;
end
% h=h_OBJ;





