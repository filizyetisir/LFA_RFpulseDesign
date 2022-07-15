


function [f df]=f0_LFA_MZ_noB0(xrf,adj,spokes,system,opt)
% function [f df]=f0_LFA_MZ_noB0(xrf,adj,spokes,system,opt)
% objective function for design of a rotation pulse (as opposed to a refocusing pulse)
% assuming the initial state Mxy(0)=0 and Mz(0)=1. Assume no B0!!

gamma=42.576*1e6;

ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);

if nargout>1
    as_gblip=ones(adj.nnonzeropixels,spokes.nspokes);
    bs_gblip=zeros(adj.nnonzeropixels,spokes.nspokes);
    
    a_forw=ones(adj.nnonzeropixels,spokes.nspokes);
    b_forw=zeros(adj.nnonzeropixels,spokes.nspokes);
    
    da_dre=cell(spokes.nspokes,1);
    da_dim=cell(spokes.nspokes,1);    
    db_dre=cell(spokes.nspokes,1);
    db_dim=cell(spokes.nspokes,1);
end

a=ones(adj.nnonzeropixels,1);
b=zeros(adj.nnonzeropixels,1);

% FORWARD BLOCH SIMULATION
for i=1:spokes.nspokes
    
    sp_ind=(i-1)*system.ncoils+1:i*system.ncoils;    
    btot=adj.b1maps_*rf(sp_ind)*spokes.sumsinc;  % total B1+ at every location
    exp_btot_pha=exp(1j*angle(btot));
    
    phi=-2*pi*gamma*system.deltat*abs(btot);  % see Pauly et al.
    cos_phi_over_2=cos(phi/2);  % no need to recompute these several times...
    sin_phi_over_2=sin(phi/2);
    
    % CK parameters incorporating gradient blips (a gradien blip is a a pure alpha with beta=0)
    a2=cos_phi_over_2 .* spokes.q_gblips(:,i);
    b2=-1j*exp_btot_pha.*sin_phi_over_2 .* conj(spokes.q_gblips(:,i));
    
    [a b]=multiply_2spinors(a,b,a2,b2);
    
    if nargout>1
        % CK parameters 
        as_gblip(:,i)=a2;
        bs_gblip(:,i)=b2;
        
        % forward Q products (1-step late)
        if i<spokes.nspokes
            a_forw(:,i+1)=a;
            b_forw(:,i+1)=b;
        end
        
        % Jacobian terms (@ each spatial position)
        cos_phi_over_2_big=repmat(cos_phi_over_2,[1 system.ncoils]);
        sin_phi_over_2_big=repmat(sin_phi_over_2,[1 system.ncoils]);
        
        tmp=adj.b1maps_.*repmat( conj(exp_btot_pha),[1 system.ncoils] );
        dphi_dre=-2*pi*gamma*(system.deltat)*abs(spokes.sumsinc)*real(tmp);
        dphi_dim=2*pi*gamma*(system.deltat)*abs(spokes.sumsinc)*imag(tmp);
        
        da_dre{i}=-0.5*sin_phi_over_2_big.*dphi_dre;
        da_dim{i}=-0.5*sin_phi_over_2_big.*dphi_dim;
        
        exp_btot_pha_big=repmat(exp_btot_pha,[1 system.ncoils]);
        btot_big=repmat(abs(btot),[1 system.ncoils]);
        dexppha_dre=( spokes.sumsinc*adj.b1maps_ - exp_btot_pha_big.*real( spokes.sumsinc*adj.b1maps_.*conj(exp_btot_pha_big) ) )./btot_big;  % derivative of the term exp(1j*angle(btot))
        dexppha_dim=( 1j*spokes.sumsinc*adj.b1maps_ - exp_btot_pha_big.*real( 1j*spokes.sumsinc*adj.b1maps_.*conj(exp_btot_pha_big) ) )./btot_big;
        dexppha_dre(abs(btot_big)==0)=0;
        dexppha_dim(abs(btot_big)==0)=0;
        
        db_dre{i}=-1j*sin_phi_over_2_big.*dexppha_dre - 0.5*1j*cos_phi_over_2_big.*exp_btot_pha_big.*dphi_dre;
        db_dim{i}=-1j*sin_phi_over_2_big.*dexppha_dim - 0.5*1j*cos_phi_over_2_big.*exp_btot_pha_big.*dphi_dim;
        
        % incorporate gradient blips in Jacobian terms
        q_gblip_big=repmat(spokes.q_gblips(:,i),[1 system.ncoils]);        
        da_dre{i}=da_dre{i} .* q_gblip_big;
        da_dim{i}=da_dim{i} .* q_gblip_big;        
        
        db_dre{i}=db_dre{i} .* conj(q_gblip_big);
        db_dim{i}=db_dim{i} .* conj(q_gblip_big);
    end
    
end


if opt.refoc_pulse==0  % excitation pulse
    if opt.do_mls==1  % MLS
        opt_metric=abs(a).^2 - abs(b).^2;  % mz=|alpha|^2-|beta|^2
    else  % LS
        opt_metric=abs(a-b).^2;  % |alpha-beta|^2
    end
else  % refocusing pulse
    if opt.do_mls==1  % MLS
        opt_metric=abs(b).^2;  % |beta|^2
    else  % LS
        opt_metric=b.^2;  % beta^2
    end
end

y=opt_metric - adj.target_LFA;
f=y'*y;


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
    
    for i=spokes.nspokes:-1:1        
        
        a_forw_repmat=repmat( a_forw(:,i),[1 system.ncoils] );
        b_forw_repmat=repmat( b_forw(:,i),[1 system.ncoils] );
        
        a_repmat=repmat( a,[1 system.ncoils] );
        b_repmat=repmat( b,[1 system.ncoils] );
        
        % update derivatives of total Q matrices (sandwitch product)
        [tmp3 tmp4]=multiply_3spinors(a_forw_repmat,b_forw_repmat,da_dre{i},db_dre{i},a_repmat,b_repmat);
        datot_dre(:,:,i)=datot_dre(:,:,i) + tmp3;
        dbtot_dre(:,:,i)=dbtot_dre(:,:,i) + tmp4;
        
        [tmp3 tmp4]=multiply_3spinors(a_forw_repmat,b_forw_repmat,da_dim{i},db_dim{i},a_repmat,b_repmat);
        datot_dim(:,:,i)=datot_dim(:,:,i) + tmp3;
        dbtot_dim(:,:,i)=dbtot_dim(:,:,i) + tmp4;
        
        % keep track of the backward Q matrices product
        [a b]=multiply_2spinors( as_gblip(:,i),bs_gblip(:,i),a,b );        
    end
    
    datot_dre=reshape( datot_dre,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    datot_dim=reshape( datot_dim,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    dbtot_dre=reshape( dbtot_dre,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    dbtot_dim=reshape( dbtot_dim,adj.nnonzeropixels,system.ncoils*spokes.nspokes );
    
    if opt.refoc_pulse==0 && opt.do_mls==0
        damb2_dre=2.0*( datot_dre-dbtot_dre ).*repmat(conj(a-b),[1 system.ncoils*spokes.nspokes]);
        damb2_dim=2.0*( datot_dim-dbtot_dim ).*repmat(conj(a-b),[1 system.ncoils*spokes.nspokes]);
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
            tmp1=real(damb2_dre);
            tmp2=real(damb2_dim);
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

















