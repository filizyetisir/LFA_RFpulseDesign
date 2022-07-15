


function [f df]=f0_LFA_mexCudaInterf(xrf,adj,spokes,system,opt)
% function [f df]=f0_LFA_MZ_mexInterf(xrf,adj,spokes,system,opt)

global gpuStruct;

gamma=42.576*1e6;

ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);
rf_exp=zeros(spokes.ntimes,system.ncoils);
for i=1:spokes.nspokes
    ind=(i-1)*system.ncoils + 1 : i*system.ncoils;
    rf_exp( spokes.sinc_time_and_spoke_to_time(:,i),: )=spokes.sinc * rf(ind).';
end

comp_grad=0;
if nargout==2
    comp_grad=1;
    dftot=zeros(2*ncunk,1);
end


if isempty(gpuStruct)  % no GPU
    
    b1s=adj.b1maps_ * rf_exp.';
    grad=adj.nonzeropixels(:,4:6) * spokes.g.' - repmat( adj.b0map_/gamma,[1 spokes.ntimes] );    
    
else  % GPU    
    
    b1s=gather( gpuArray(adj.b1maps_) * gpuArray(rf_exp.') );
    grad=gather( gpuArray(adj.nonzeropixels(:,4:6)) * gpuArray(spokes.g.') ) - repmat( adj.b0map_/gamma,[1 spokes.ntimes] );        
    
end


ftot=0;

for i=1:opt.nfreqs    
    
    grad2=grad - opt.freqs(i)/gamma;    
    
    if isempty(gpuStruct)  % no GPU --> use MEX function
        
        inds=(i-1)*adj.nnonzeropixels + 1 : i*adj.nnonzeropixels;
        
        if(opt.multiCPU == 1)
            % If you want to parallelize over all CPU processors
            if(opt.useJ == 1)
                comp_grad = 1;
                [f df]=f0_LFA_mex(b1s,grad2,spokes.time_to_sinc_time,spokes.time_to_spoke, ...
                    spokes.sinc,adj.b1maps_,complex(real(adj.target_LFA(inds)),imag(adj.target_LFA(inds))),spokes.nspokes,system.deltat,comp_grad,opt.do_mls,opt.refoc_pulse);
            else
                comp_grad = 0;
                [f]=f0_LFA_mex(b1s,grad2,spokes.time_to_sinc_time,spokes.time_to_spoke, ...
                    spokes.sinc,adj.b1maps_,complex(real(adj.target_LFA(inds)),imag(adj.target_LFA(inds))),spokes.nspokes,system.deltat,comp_grad,opt.do_mls,opt.refoc_pulse);
            end
        else
            % If you want to use only single CPU processor
            if(opt.useJ == 1)
                [f df]=f0_LFA(xrf,adj,spokes,system,opt);
            else
                [f]=f0_LFA(xrf,adj,spokes,system,opt);
            end
        end
        
    else  % GPU --> use GPU kernel
        
        % fast forward and backward Bloch sims using GPU kernel
         [gpuStruct.a_re gpuStruct.a_im gpuStruct.b_re gpuStruct.b_im ...
          gpuStruct.datot_dre_re gpuStruct.datot_dre_im gpuStruct.datot_dim_re gpuStruct.datot_dim_im ...
          gpuStruct.dbtot_dre_re gpuStruct.dbtot_dre_im gpuStruct.dbtot_dim_re gpuStruct.dbtot_dim_im ]=feval( gpuStruct.f0_LFA_gpuKernel, ...
             gpuStruct.a_re, gpuStruct.a_im, gpuStruct.b_re, gpuStruct.b_im, ...				  
             gpuStruct.datot_dre_re, gpuStruct.datot_dre_im, gpuStruct.datot_dim_re, gpuStruct.datot_dim_im , ...
             gpuStruct.dbtot_dre_re, gpuStruct.dbtot_dre_im, gpuStruct.dbtot_dim_re, gpuStruct.dbtot_dim_im , ...
             gpuStruct.a_forw_re ,gpuStruct.a_forw_im ,gpuStruct.b_forw_re ,gpuStruct.b_forw_im , ...
             real(b1s),imag(b1s),grad2, ...
             spokes.time_to_sinc_time,spokes.time_to_spoke,spokes.sinc, ...
             real(adj.b1maps_),imag(adj.b1maps_), ...
             adj.nnonzeropixels,spokes.ntimes,system.ncoils,spokes.nspokes,system.deltat,comp_grad );

         a=gpuArray( gpuStruct.a_re + 1j*gpuStruct.a_im );
         b=gpuArray( gpuStruct.b_re + 1j*gpuStruct.b_im );
         
         % objective function
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
         
         inds=(i-1)*adj.nnonzeropixels + 1 : i*adj.nnonzeropixels;
         y=opt_metric - adj.target_LFA(inds);
         f=real( gather( y'*y ) );
         
         % gradient
         if comp_grad==1
             
             datot_dre=gpuArray( reshape( gpuStruct.datot_dre_re + 1j*gpuStruct.datot_dre_im,adj.nnonzeropixels,system.ncoils*spokes.nspokes ) );
             datot_dim=gpuArray( reshape( gpuStruct.datot_dim_re + 1j*gpuStruct.datot_dim_im,adj.nnonzeropixels,system.ncoils*spokes.nspokes ) );    
             dbtot_dre=gpuArray( reshape( gpuStruct.dbtot_dre_re + 1j*gpuStruct.dbtot_dre_im,adj.nnonzeropixels,system.ncoils*spokes.nspokes ) );
             dbtot_dim=gpuArray( reshape( gpuStruct.dbtot_dim_re + 1j*gpuStruct.dbtot_dim_im,adj.nnonzeropixels,system.ncoils*spokes.nspokes ) );
             
             if opt.refoc_pulse==0 && opt.do_mls==0
                dab1_dre=2.0*conj(datot_dre) .*repmat(b,[1 system.ncoils*spokes.nspokes]);
                dab1_dim=2.0*conj(datot_dim) .*repmat(b,[1 system.ncoils*spokes.nspokes]);

                dab2_dre=2.0*dbtot_dre .*repmat(conj(a),[1 system.ncoils*spokes.nspokes]);
                dab2_dim=2.0*dbtot_dim .*repmat(conj(a),[1 system.ncoils*spokes.nspokes]);
             elseif opt.refoc_pulse==1 && opt.do_mls==0
                 da2_dre=2.0*datot_dre.*repmat(a,[1 system.ncoils*spokes.nspokes]);
                 da2_dim=2.0*datot_dim.*repmat(a,[1 system.ncoils*spokes.nspokes]);
                 
                 db2_dre=2.0*dbtot_dre.*repmat(b,[1 system.ncoils*spokes.nspokes]);
                 db2_dim=2.0*dbtot_dim.*repmat(b,[1 system.ncoils*spokes.nspokes]);
             else
                 da2_dre=2.0*datot_dre.*repmat(conj(a),[1 system.ncoils*spokes.nspokes]);
                 da2_dim=2.0*datot_dim.*repmat(conj(a),[1 system.ncoils*spokes.nspokes]);
                 
                 db2_dre=2.0*dbtot_dre.*repmat(conj(b),[1 system.ncoils*spokes.nspokes]);
                 db2_dim=2.0*dbtot_dim.*repmat(conj(b),[1 system.ncoils*spokes.nspokes]);
             end
             
             if opt.refoc_pulse==0  % excitation pulse
                 if opt.do_mls==1  % MLS
                     tmp1=real(da2_dre) - real(db2_dre);
                     tmp2=real(da2_dim) - real(db2_dim);
                 else  % LS
                    tmp1=dab1_dre + dab2_dre;
                    tmp2=dab1_dim + dab2_dim;
                 end
             else  % refocusing pulse
                 if opt.do_mls==1  % MLS
                     tmp1=real(db2_dre);
                     tmp2=real(db2_dim);
                 else  % LS
                     tmp1=db2_dre;
                     tmp2=db2_dim;
                 end
             end
             
             df=2.0*real( gather( y'*[tmp1 tmp2] ) )';
             
         end         
    end
    
    % total objective function is the sum of obj fct at optimized
    % frequencies
    ftot=ftot + f*opt.freq_weights(i);
    
    if comp_grad==1
        dftot=dftot + df*opt.freq_weights(i);
    end
    
end

f = ftot;
if(comp_grad==1)
df = dftot;
end

































    




