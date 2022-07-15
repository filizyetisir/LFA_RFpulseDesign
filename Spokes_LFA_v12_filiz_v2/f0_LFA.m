


function [f df]=f0_LFA(xrf,adj,spokes,system,opt)
% function [f df]=f0_LFA(xrf,adj,spokes,system,opt)
% objective function (with Jacobian is requested) for design of excitation and refocusing pulses. LS and
% MLS designs are supported

gamma=42.576*1e6;

ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);
rf_exp=zeros(spokes.ntimes,system.ncoils);
for i=1:spokes.nspokes
    ind=(i-1)*system.ncoils + 1 : i*system.ncoils;
    rf_exp( spokes.sinc_time_and_spoke_to_time(:,i),: )=spokes.sinc * rf(ind).';
end

if nargout>1
    % if gradient is requested, store all forward
    % Q matrix products
    a_forw=ones(adj.nnonzeropixels,spokes.ntimes); 
    b_forw=zeros(adj.nnonzeropixels,spokes.ntimes);
    
    % actual Q matrices
    as=zeros(adj.nnonzeropixels,spokes.ntimes); 
    bs=zeros(adj.nnonzeropixels,spokes.ntimes);
    
    % derivatives of Q matrices
    da_dre=zeros(spokes.ntimes,adj.nnonzeropixels,system.ncoils);
    da_dim=zeros(spokes.ntimes,adj.nnonzeropixels,system.ncoils);    
    db_dre=zeros(spokes.ntimes,adj.nnonzeropixels,system.ncoils);
    db_dim=zeros(spokes.ntimes,adj.nnonzeropixels,system.ncoils);
end

a=ones(adj.nnonzeropixels,1);
b=zeros(adj.nnonzeropixels,1);

% FORWARD BLOCH SIMULATION
for time=1:spokes.ntimes
    
    b1s=( rf_exp(time,:)*adj.b1maps_.' ).';
    grad=adj.nonzeropixels(:,4:6) * spokes.g(time,:).' - adj.b0map_/gamma;
    norm=sqrt( abs(b1s).^2+grad.^2 );
    phi=-2*pi*gamma*system.deltat*norm;
    
    n=[real(b1s) imag(b1s) grad];
    n=n./repmat(norm,1,3);
    
    a2=cos(phi/2.0) - 1j*(n(:,3)).*sin(phi/2.0);
    b2=-1j*( n(:,1)+1j*n(:,2) ).*sin(phi/2.0);
       
    ind2=find(norm==0 | isnan(norm));
    a2(ind2)=1.0;
    b2(ind2)=0.0;
    
    % Q matrix product
    [a b]=multiply_2spinors(a,b,a2,b2);

    if nargout>1
        
        % store Q matrix
        as(:,time)=a2;
        bs(:,time)=b2;
        
        % store one time late forward Q matrix product
        if time<spokes.ntimes
            a_forw(:,time+1)=a;
            b_forw(:,time+1)=b;
        end
               
        sinctime=spokes.time_to_sinc_time(time);
        if sinctime>0
            % derivative of phi wrt spokes amplitudes
            dphi=-2*pi*gamma*system.deltat .* repmat(conj(b1s)./norm,[1 system.ncoils]) .* spokes.sinc(sinctime) .* adj.b1maps_;
            dphi_dre=real(dphi);
            dphi_dim=-imag(dphi);
            
            % derivative of n (rotation axis) wrt spokes amplitudes
            tmp=-repmat(conj(b1s)./norm.^3,[1 system.ncoils]) .* spokes.sinc(sinctime) .* adj.b1maps_ ;
            tmp2=repmat(b1s,[1 system.ncoils]);
            tmp3=repmat( spokes.sinc(sinctime)./norm,[1 system.ncoils] );
            
            dnx_dre=real(tmp).*real(tmp2)  + tmp3.*real(adj.b1maps_);
            dnx_dim=-imag(tmp).*real(tmp2) - tmp3.*imag(adj.b1maps_);
            
            dny_dre=real(tmp).*imag(tmp2) + tmp3.*imag(adj.b1maps_);
            dny_dim=-imag(tmp).*imag(tmp2) + tmp3.*real(adj.b1maps_);

            dnz=tmp .* repmat(grad,[1 system.ncoils]);
            dnz_dre=real(dnz);
            dnz_dim=-imag(dnz);            
            
            % derivative of Q
            da_dre(time,:,:)=dphi_dre/2.*repmat( -sin(phi/2) - 1j*n(:,3).*cos(phi/2),[1 system.ncoils] ) - 1j*repmat(sin(phi/2),[1 system.ncoils]).*dnz_dre;
            da_dim(time,:,:)=dphi_dim/2.*repmat( -sin(phi/2) - 1j*n(:,3).*cos(phi/2),[1 system.ncoils] ) - 1j*repmat(sin(phi/2),[1 system.ncoils]).*dnz_dim;
            
            db_dre(time,:,:)=-1j*( dnx_dre + 1j*dny_dre ) .* repmat(sin(phi/2),[1 system.ncoils]) - 1j/2*repmat( ( n(:,1)+1j*n(:,2) ).*cos(phi/2),[1 system.ncoils] ) .* dphi_dre;        
            db_dim(time,:,:)=-1j*( dnx_dim + 1j*dny_dim ) .* repmat(sin(phi/2),[1 system.ncoils]) - 1j/2*repmat( ( n(:,1)+1j*n(:,2) ).*cos(phi/2),[1 system.ncoils] ) .* dphi_dim;
            
            % remove problematic indices from derivatives
            norm2=repmat(norm,[1 system.ncoils]);
            ind3=find( norm2==0 | isnan(norm2) );
            
            da_dre(time,ind3)=0;
            da_dim(time,ind3)=0;            
            db_dre(time,ind3)=0;
            db_dim(time,ind3)=0;        
        else
            da_dre(time,:)=0;
            da_dim(time,:)=0;            
            db_dre(time,:)=0;
            db_dim(time,:)=0;
        end
        
    end
    
end

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

y=opt_metric - adj.target_LFA;
f=real( y'*y );  % remove small imagingary values due to numerical errors


% BACKWARD BLOCH SIMULATION
if nargout>1
    
    da_dre=permute( da_dre,[2 3 1] );
    da_dim=permute( da_dim,[2 3 1] );

    db_dre=permute( db_dre,[2 3 1] );
    db_dim=permute( db_dim,[2 3 1] );

    a=ones(adj.nnonzeropixels,1);
    b=zeros(adj.nnonzeropixels,1);
    
    datot_dre=zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes);
    datot_dim=zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes);
    dbtot_dre=zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes);
    dbtot_dim=zeros(adj.nnonzeropixels,system.ncoils,spokes.nspokes);
    
    for time=spokes.ntimes:-1:1        
        
        sp_num=spokes.time_to_spoke(time);
        if sp_num>0
            a_forw_repmat=repmat( a_forw(:,time),[1 system.ncoils] );
            b_forw_repmat=repmat( b_forw(:,time),[1 system.ncoils] );
            
            a_repmat=repmat( a,[1 system.ncoils] );
            b_repmat=repmat( b,[1 system.ncoils] );
            
            % update derivatives of total Q matrices (sandwitch product)
            [tmp3 tmp4]=multiply_3spinors(a_forw_repmat,b_forw_repmat,da_dre(:,:,time),db_dre(:,:,time),a_repmat,b_repmat);
            datot_dre(:,:,sp_num)=datot_dre(:,:,sp_num) + tmp3;
            dbtot_dre(:,:,sp_num)=dbtot_dre(:,:,sp_num) + tmp4;
            
            [tmp3 tmp4]=multiply_3spinors(a_forw_repmat,b_forw_repmat,da_dim(:,:,time),db_dim(:,:,time),a_repmat,b_repmat);
            datot_dim(:,:,sp_num)=datot_dim(:,:,sp_num) + tmp3;
            dbtot_dim(:,:,sp_num)=dbtot_dim(:,:,sp_num) + tmp4;
        end
        
        % keep track of the backward Q matrices product
        [a b]=multiply_2spinors( as(:,time),bs(:,time),a,b );        
    end
    
    datot_dre=reshape( datot_dre,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    datot_dim=reshape( datot_dim,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    dbtot_dre=reshape( dbtot_dre,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    dbtot_dim=reshape( dbtot_dim,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    
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
    
    df=2.0*real( y'*[tmp1 tmp2] )';
end
















