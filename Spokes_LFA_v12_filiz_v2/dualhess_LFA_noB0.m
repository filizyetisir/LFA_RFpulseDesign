


function h=dualhess_LFA_noB0(xrf,lambda,adj,system,spokes,opt,sar,solver)
% function h=dualhess_LFA_noB0(xrf,lambda,adj,system,spokes,opt,sar,solver)
% Works for any number of spokes. 

gamma=42.576*1e6;

ncunk=size(xrf,1)/2;
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);


% FORWARD BLOCH SIMULATION
as_gblip=ones(adj.nnonzeropixels,spokes.nspokes);
bs_gblip=zeros(adj.nnonzeropixels,spokes.nspokes);

a_forw=ones(adj.nnonzeropixels,spokes.nspokes);
b_forw=zeros(adj.nnonzeropixels,spokes.nspokes);

h_a_1=cell(spokes.nspokes,1);
h_a_2=cell(spokes.nspokes,1);
h_a_3=cell(spokes.nspokes,1);

da_dre=cell(spokes.nspokes,1);
da_dim=cell(spokes.nspokes,1);

h_b_1=cell(spokes.nspokes,1);
h_b_2=cell(spokes.nspokes,1);
h_b_3=cell(spokes.nspokes,1);

db_dre=cell(spokes.nspokes,1);
db_dim=cell(spokes.nspokes,1);

for i=1:spokes.nspokes
    
    sp_ind=(i-1)*system.ncoils+1:i*system.ncoils;    
    btot=adj.b1maps_*rf(sp_ind)*spokes.sumsinc;  % total B1+ at every location
    exp_btot_pha=exp(1j*angle(btot));
    
    phi=-2.0*pi*gamma*system.deltat*abs(btot);  % see Pauly et al.
    cos_phi_over_2=cos(phi/2);  % no need to recompute these several times...
    sin_phi_over_2=sin(phi/2);    
    cos_phi_over_2_big=repmat(cos_phi_over_2,[1 system.ncoils]);
    sin_phi_over_2_big=repmat(sin_phi_over_2,[1 system.ncoils]);
        
    % CK parameters
    a=cos_phi_over_2;
    b=-1j*exp_btot_pha.*sin_phi_over_2;
    
    % CK parameters incorporating gradient blips (a gradien blip is a a pure alpha with beta=0)
    as_gblip(:,i)=a.*spokes.q_gblips(:,i);
    bs_gblip(:,i)=b.*conj(spokes.q_gblips(:,i));
    
    % forward Q products (1-step late)
    if i<spokes.nspokes
        [a_forw(:,i+1) b_forw(:,i+1)]=multiply_2spinors(a_forw(:,i),b_forw(:,i),as_gblip(:,i),bs_gblip(:,i));    
    end

    % Jacobian terms (@ each spatial position)
    tmp=adj.b1maps_.*repmat( conj(exp_btot_pha),[1 system.ncoils] );
    dphi_dre=-2.0*pi*gamma*(system.deltat)*abs(spokes.sumsinc)*real(tmp);
    dphi_dim=2.0*pi*gamma*(system.deltat)*abs(spokes.sumsinc)*imag(tmp);
        
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
    
    % Hessian of phi (@ each spatial position)
    norm=-2.0*pi*gamma*(system.deltat)*(spokes.sumsinc)^2./abs(btot);
    norm(abs(btot)==0)=0;
    h_phi_1=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_phi_2=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_phi_3=zeros(adj.nnonzeropixels,system.ncoils^2);
    tmp2=adj.b1maps_.*conj(exp_btot_pha_big);
    for j=1:system.ncoils
        for k=1:system.ncoils
            ind=(j-1)*system.ncoils + k;
            tmp=adj.b1maps_(:,j).*conj(adj.b1maps_(:,k));
            if k>=j
                h_phi_1(:,ind)=norm.*( real(tmp) - real(tmp2(:,j)).*real(tmp2(:,k)) );            
                h_phi_3(:,ind)=norm.*( real(tmp) - imag(tmp2(:,j)).*imag(tmp2(:,k)) );                                
                if k~=j
                    ind2=(k-1)*system.ncoils + j;
                    h_phi_1(:,ind2)=h_phi_1(:,ind);  % symmetric term
                    h_phi_3(:,ind2)=h_phi_3(:,ind);  % symmetric term
                end
            end
            h_phi_2(:,ind)=norm.*( -imag(tmp) + imag(tmp2(:,j)).*real(tmp2(:,k)) );  % not symmetric
        end
    end
    
    
    % Hessian of a (@ each spatial location)
    h_a_1{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_a_2{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_a_3{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    for j=1:system.ncoils
        for k=1:system.ncoils
            ind=(j-1)*system.ncoils + k;
            if k>=j
                h_a_1{i}(:,ind)=-0.25*cos_phi_over_2.*dphi_dre(:,j).*dphi_dre(:,k) - 0.5*sin_phi_over_2.*h_phi_1(:,ind);  % hessian term #1 (dre_dre)
                h_a_3{i}(:,ind)=-0.25*cos_phi_over_2.*dphi_dim(:,j).*dphi_dim(:,k) - 0.5*sin_phi_over_2.*h_phi_3(:,ind);  % hessian term #3 (dim_dim)
                if k~=j
                    ind2=(k-1)*system.ncoils + j;
                    h_a_1{i}(:,ind2)=h_a_1{i}(:,ind);  % symmetric term
                    h_a_3{i}(:,ind2)=h_a_3{i}(:,ind);  % symmetric term
                end
            end
            h_a_2{i}(:,ind)=-0.25*cos_phi_over_2.*dphi_dim(:,j).*dphi_dre(:,k) - 0.5*sin_phi_over_2.*h_phi_2(:,ind);  % hessian term #2 (dre_dim)
        end
    end
    
    % Hessian of the term exp(1j*angle(btot)) (@ each spatial position)
    h_exppha_1=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_exppha_2=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_exppha_3=zeros(adj.nnonzeropixels,system.ncoils^2);
    ind_null=find(abs(btot)==0);
    V=abs(btot);
    for j=1:system.ncoils
        for k=1:system.ncoils
            ind=(j-1)*system.ncoils + k;
            if k>=j
                % term #1
                U=spokes.sumsinc*( adj.b1maps_(:,j) - exp_btot_pha.*real( adj.b1maps_(:,j).*conj(exp_btot_pha) ) );                
                dU=-dexppha_dre(:,k).*real( spokes.sumsinc*adj.b1maps_(:,j).*conj(exp_btot_pha) ) - exp_btot_pha.*real( spokes.sumsinc*adj.b1maps_(:,j).*conj(dexppha_dre(:,k) ) );
                dV=real( spokes.sumsinc*adj.b1maps_(:,k).*conj(exp_btot_pha) );
                h_exppha_1(:,ind)=( dU.*V-dV.*U )./V.^2;
                h_exppha_1(ind_null,ind)=0;

                % term #3
                U=spokes.sumsinc*( 1j*adj.b1maps_(:,j) - exp_btot_pha.*real( 1j*adj.b1maps_(:,j).*conj(exp_btot_pha) ) );
                dU=dexppha_dim(:,k).*imag( spokes.sumsinc*adj.b1maps_(:,j).*conj(exp_btot_pha) ) + exp_btot_pha.*imag( spokes.sumsinc*adj.b1maps_(:,j).*conj(dexppha_dim(:,k) ) );
                dV=-imag( spokes.sumsinc*adj.b1maps_(:,k).*conj(exp_btot_pha) );
                h_exppha_3(:,ind)=( dU.*V-dV.*U )./V.^2;
                h_exppha_3(ind_null,ind)=0;

                if k~=j
                    ind2=(k-1)*system.ncoils + j;
                    h_exppha_1(:,ind2)=h_exppha_1(:,ind);  % symmetric term
                    h_exppha_3(:,ind2)=h_exppha_3(:,ind);  % symmetric term
                end
            end
            
            % term #2
            U=spokes.sumsinc*( 1j*adj.b1maps_(:,j) - exp_btot_pha.*real( 1j*adj.b1maps_(:,j).*conj(exp_btot_pha) ) );
            dU=dexppha_dre(:,k).*imag( spokes.sumsinc*adj.b1maps_(:,j).*conj(exp_btot_pha) ) + exp_btot_pha.*imag( spokes.sumsinc*adj.b1maps_(:,j).*conj(dexppha_dre(:,k) ) );
            dV=real( spokes.sumsinc*adj.b1maps_(:,k).*conj(exp_btot_pha) );
            h_exppha_2(:,ind)=( dU.*V-dV.*U )./V.^2;
            h_exppha_2(ind_null,ind)=0;            
        end
    end
    
    % Hessian of b (@ each spatial position)
    h_b_1{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_b_2{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    h_b_3{i}=zeros(adj.nnonzeropixels,system.ncoils^2);
    for j=1:system.ncoils
        for k=1:system.ncoils
            ind=(j-1)*system.ncoils + k;
            if k>=j
                % term #1
                tmp1=-0.5*1j*cos_phi_over_2.*( dphi_dre(:,k).*dexppha_dre(:,j)+dphi_dre(:,j).*dexppha_dre(:,k) );
                tmp2=-1j*sin_phi_over_2.*h_exppha_1(:,ind);
                tmp3=0.25*1j*sin_phi_over_2.*exp_btot_pha.*dphi_dre(:,k).*dphi_dre(:,j);
                tmp4=-0.5*1j*cos_phi_over_2.*exp_btot_pha.*h_phi_1(:,ind);
                h_b_1{i}(:,ind)=tmp1+tmp2+tmp3+tmp4;

                % term #3
                tmp1=-0.5*1j*cos_phi_over_2.*( dphi_dim(:,k).*dexppha_dim(:,j)+dphi_dim(:,j).*dexppha_dim(:,k) );
                tmp2=-1j*sin_phi_over_2.*h_exppha_3(:,ind);
                tmp3=0.25*1j*sin_phi_over_2.*exp_btot_pha.*dphi_dim(:,k).*dphi_dim(:,j);
                tmp4=-0.5*1j*cos_phi_over_2.*exp_btot_pha.*h_phi_3(:,ind);
                h_b_3{i}(:,ind)=tmp1+tmp2+tmp3+tmp4;

                if k~=j
                    ind2=(k-1)*system.ncoils + j;
                    h_b_1{i}(:,ind2)=h_b_1{i}(:,ind);  % symmetric term
                    h_b_3{i}(:,ind2)=h_b_3{i}(:,ind);  % symmetric term
                end
            end
            
            % term #2
            tmp1=-0.5*1j*cos_phi_over_2.*( dphi_dre(:,k).*dexppha_dim(:,j)+dphi_dim(:,j).*dexppha_dre(:,k) );
            tmp2=-1j*sin_phi_over_2.*h_exppha_2(:,ind);
            tmp3=0.25*1j*sin_phi_over_2.*exp_btot_pha.*dphi_dre(:,k).*dphi_dim(:,j);
            tmp4=-0.5*1j*cos_phi_over_2.*exp_btot_pha.*h_phi_2(:,ind);
            h_b_2{i}(:,ind)=tmp1+tmp2+tmp3+tmp4;
        end
    end
    
    % incorporate gradient blips in Hessian and Jacobian terms
    q_gblip_big=repmat(spokes.q_gblips(:,i),[1 system.ncoils^2]);
    h_a_1{i}=h_a_1{i} .* q_gblip_big;
    h_a_2{i}=h_a_2{i} .* q_gblip_big;
    h_a_3{i}=h_a_3{i} .* q_gblip_big;
    
    da_dre{i}=da_dre{i} .* q_gblip_big(:,1:system.ncoils);
    da_dim{i}=da_dim{i} .* q_gblip_big(:,1:system.ncoils);
    
    h_b_1{i}=h_b_1{i} .* conj(q_gblip_big);
    h_b_2{i}=h_b_2{i} .* conj(q_gblip_big);
    h_b_3{i}=h_b_3{i} .* conj(q_gblip_big);
        
    db_dre{i}=db_dre{i} .* conj(q_gblip_big(:,1:system.ncoils));
    db_dim{i}=db_dim{i} .* conj(q_gblip_big(:,1:system.ncoils));    
end



% BACKWARD BLOCH SIMULATION
a=ones(adj.nnonzeropixels,1);  
b=zeros(adj.nnonzeropixels,1);

ncunk=spokes.nspokes*system.ncoils;

ha=zeros(adj.nnonzeropixels,2.0*ncunk,2.0*ncunk);    
hb=zeros(adj.nnonzeropixels,2.0*ncunk,2.0*ncunk);

ja_re=zeros(adj.nnonzeropixels,ncunk);
ja_im=zeros(adj.nnonzeropixels,ncunk);

jb_re=zeros(adj.nnonzeropixels,ncunk);
jb_im=zeros(adj.nnonzeropixels,ncunk);

for i=spokes.nspokes:-1:1

    ind=(i-1)*system.ncoils+1 : i*system.ncoils;
        
    a_big=repmat(a,[1 system.ncoils^2]);
    b_big=repmat(b,[1 system.ncoils^2]);
    a_forw_big=repmat(a_forw(:,i),[1 system.ncoils^2]);
    b_forw_big=repmat(b_forw(:,i),[1 system.ncoils^2]);
    
    % Jacobian (sandwitch product of backward Q product, Jacobian of current Q matrix and forward Q products)
    [ja_re(:,ind) jb_re(:,ind)]=multiply_3spinors(a_forw_big(:,1:system.ncoils),b_forw_big(:,1:system.ncoils),da_dre{i},db_dre{i},a_big(:,1:system.ncoils),b_big(:,1:system.ncoils));
    [ja_im(:,ind) jb_im(:,ind)]=multiply_3spinors(a_forw_big(:,1:system.ncoils),b_forw_big(:,1:system.ncoils),da_dim{i},db_dim{i},a_big(:,1:system.ncoils),b_big(:,1:system.ncoils));

    % Hessian term (derivatives taken within a single spoke)
    [ha1 hb1]=multiply_3spinors(a_forw_big,b_forw_big,h_a_1{i},h_b_1{i},a_big,b_big);
    [ha2 hb2]=multiply_3spinors(a_forw_big,b_forw_big,h_a_2{i},h_b_2{i},a_big,b_big);
    [ha3 hb3]=multiply_3spinors(a_forw_big,b_forw_big,h_a_3{i},h_b_3{i},a_big,b_big);
    
    ha(:,ind,ind)=reshape(ha1,adj.nnonzeropixels,system.ncoils,system.ncoils);
    ha(:,ind,ind+ncunk)=reshape(ha2,adj.nnonzeropixels,system.ncoils,system.ncoils);
    ha(:,ind+ncunk,ind)=permute( ha(:,ind,ind+ncunk),[1 3 2] );
    ha(:,ind+ncunk,ind+ncunk)=reshape(ha3,adj.nnonzeropixels,system.ncoils,system.ncoils);

    hb(:,ind,ind)=reshape(hb1,adj.nnonzeropixels,system.ncoils,system.ncoils);
    hb(:,ind,ind+ncunk)=reshape(hb2,adj.nnonzeropixels,system.ncoils,system.ncoils);
    hb(:,ind+ncunk,ind)=permute( hb(:,ind,ind+ncunk),[1 3 2] );
    hb(:,ind+ncunk,ind+ncunk)=reshape(hb3,adj.nnonzeropixels,system.ncoils,system.ncoils);
    
    % Jacobian terms (derivatives taken for different spokes)
    if spokes.nspokes>1        
        
        % product between left-hand jacobian term and backward Q products
        [a_tmp1_re b_tmp1_re]=multiply_2spinors(da_dre{i},db_dre{i},a_big(:,1:system.ncoils),b_big(:,1:system.ncoils));
        [a_tmp1_im b_tmp1_im]=multiply_2spinors(da_dim{i},db_dim{i},a_big(:,1:system.ncoils),b_big(:,1:system.ncoils));
        
        for j=i-1:-1:1            
            
            % compute the product of Q matrices sandwitched between indices i and j
            a_sandw=ones(adj.nnonzeropixels,1);
            b_sandw=zeros(adj.nnonzeropixels,1);
            for k=i-1:-1:j+1
                [a_sandw b_sandw]=multiply_2spinors(as_gblip(:,k),bs_gblip(:,k),a_sandw,b_sandw);
            end
            
            % product between sandwictch term, right-hand jacobian term and forward Q products
            a_forw_big2=repmat(a_forw(:,j),[1 system.ncoils]);
            b_forw_big2=repmat(b_forw(:,j),[1 system.ncoils]);
            [a_tmp2_re b_tmp2_re]=multiply_3spinors(a_forw_big2,b_forw_big2,da_dre{j},db_dre{j},repmat(a_sandw,[1 system.ncoils]),repmat(b_sandw,[1 system.ncoils]));
            [a_tmp2_im b_tmp2_im]=multiply_3spinors(a_forw_big2,b_forw_big2,da_dim{j},db_dim{j},repmat(a_sandw,[1 system.ncoils]),repmat(b_sandw,[1 system.ncoils]));
            
            % outter products between the columns of the TMP1 and TMP2 Q terms
            ja1=zeros(adj.nnonzeropixels,system.ncoils^2); jb1=zeros(adj.nnonzeropixels,system.ncoils^2);
            ja2=zeros(adj.nnonzeropixels,system.ncoils^2); jb2=zeros(adj.nnonzeropixels,system.ncoils^2);
            ja3=zeros(adj.nnonzeropixels,system.ncoils^2); jb3=zeros(adj.nnonzeropixels,system.ncoils^2);
            ja4=zeros(adj.nnonzeropixels,system.ncoils^2); jb4=zeros(adj.nnonzeropixels,system.ncoils^2);
            
            for i2=1:system.ncoils
                for j2=1:system.ncoils
                    ind3=(i2-1)*system.ncoils + j2;
                    [ja1(:,ind3) jb1(:,ind3)]=multiply_2spinors(a_tmp2_re(:,j2),b_tmp2_re(:,j2) , a_tmp1_re(:,i2),b_tmp1_re(:,i2));
                    [ja2(:,ind3) jb2(:,ind3)]=multiply_2spinors(a_tmp2_re(:,j2),b_tmp2_re(:,j2) , a_tmp1_im(:,i2),b_tmp1_im(:,i2));
                    [ja3(:,ind3) jb3(:,ind3)]=multiply_2spinors(a_tmp2_im(:,j2),b_tmp2_im(:,j2) , a_tmp1_re(:,i2),b_tmp1_re(:,i2));
                    [ja4(:,ind3) jb4(:,ind3)]=multiply_2spinors(a_tmp2_im(:,j2),b_tmp2_im(:,j2) , a_tmp1_im(:,i2),b_tmp1_im(:,i2));
                end
            end
            
            ind2=(j-1)*system.ncoils+1 : j*system.ncoils;

            ha(:,ind2,ind)=reshape(ja1,adj.nnonzeropixels,system.ncoils,system.ncoils);
            ha(:,ind2,ind+ncunk)=reshape(ja2,adj.nnonzeropixels,system.ncoils,system.ncoils);
            ha(:,ind2+ncunk,ind)=reshape(ja3,adj.nnonzeropixels,system.ncoils,system.ncoils);
            ha(:,ind2+ncunk,ind+ncunk)=reshape(ja4,adj.nnonzeropixels,system.ncoils,system.ncoils);
            
            hb(:,ind2,ind)=reshape(jb1,adj.nnonzeropixels,system.ncoils,system.ncoils);
            hb(:,ind2,ind+ncunk)=reshape(jb2,adj.nnonzeropixels,system.ncoils,system.ncoils);
            hb(:,ind2+ncunk,ind)=reshape(jb3,adj.nnonzeropixels,system.ncoils,system.ncoils);
            hb(:,ind2+ncunk,ind+ncunk)=reshape(jb4,adj.nnonzeropixels,system.ncoils,system.ncoils);

            % fill half of the hessian matrix by symmetry
            ha(:,ind,ind2)=permute( ha(:,ind2,ind),[1 3 2] );
            ha(:,ind+ncunk,ind2)=permute( ha(:,ind2,ind+ncunk),[1 3 2] );
            ha(:,ind,ind2+ncunk)=permute( ha(:,ind2+ncunk,ind),[1 3 2] );
            ha(:,ind+ncunk,ind2+ncunk)=permute( ha(:,ind2+ncunk,ind+ncunk),[1 3 2] );
            
            hb(:,ind,ind2)=permute( hb(:,ind2,ind),[1 3 2] );
            hb(:,ind+ncunk,ind2)=permute( hb(:,ind2,ind+ncunk),[1 3 2] );
            hb(:,ind,ind2+ncunk)=permute( hb(:,ind2+ncunk,ind),[1 3 2] );
            hb(:,ind+ncunk,ind2+ncunk)=permute( hb(:,ind2+ncunk,ind+ncunk),[1 3 2] );
        end  % for j=i-1:-1:1
    
    end  % if spokes.nspokes>1  
        
    % keep track of the backward Q products
    [a b]=multiply_2spinors(as_gblip(:,i),bs_gblip(:,i),a,b);    
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

h11=tmp1'*tmp1;
h12=tmp1'*tmp2;
h22=tmp2'*tmp2;

h1=2.0*real([ h11   h12;
              h12.' h22 ] );
     
% Hessian terms of total Hessian
ha_=reshape(ha,adj.nnonzeropixels,(2*ncunk)^2);
hb_=reshape(hb,adj.nnonzeropixels,(2*ncunk)^2);

if opt.refoc_pulse==0  % excitation pulse
    if opt.do_mls==1  % MLS
        tmp1=reshape( (y.*a)'*ha_ - (y.*b)'*hb_,2.0*ncunk,2.0*ncunk );
        tmp2=ja_re'*(repmat(y,[1 ncunk]).*ja_re) - jb_re'*(repmat(y,[1 ncunk]).*jb_re);
        tmp3=ja_im'*(repmat(y,[1 ncunk]).*ja_re) - jb_im'*(repmat(y,[1 ncunk]).*jb_re);
        tmp4=ja_im'*(repmat(y,[1 ncunk]).*ja_im) - jb_im'*(repmat(y,[1 ncunk]).*jb_im);
    else  % LS
        tmp1=reshape( (y.*conj(b))'*conj(ha_) + (y.*a)'*hb_,2.0*ncunk,2.0*ncunk );        
        tmp2=ja_re'*(repmat(conj(y),[1 ncunk]).*jb_re) + jb_re.'*(repmat(conj(y),[1 ncunk]).*conj(ja_re));
        tmp3=ja_im'*(repmat(conj(y),[1 ncunk]).*jb_re) + jb_im.'*(repmat(conj(y),[1 ncunk]).*conj(ja_re));
        tmp4=ja_im'*(repmat(conj(y),[1 ncunk]).*jb_im) + jb_im.'*(repmat(conj(y),[1 ncunk]).*conj(ja_im));        
    end
else  % refocusing pulse
    if opt.do_mls==1  % MLS
        tmp1=reshape( (y.*conj(b))'*hb_,2.0*ncunk,2.0*ncunk );
        tmp2=2.0*jb_re'*(repmat(y,[1 ncunk]).*jb_re);
        tmp3=2.0*jb_im'*(repmat(y,[1 ncunk]).*jb_re);
        tmp4=2.0*jb_im'*(repmat(y,[1 ncunk]).*jb_im);
    else  % LS
        tmp1=reshape( (y.*conj(b))'*hb_,2.0*ncunk,2.0*ncunk );
        tmp2=jb_re.'*(repmat(conj(y),[1 ncunk]).*jb_re);
        tmp3=jb_im.'*(repmat(conj(y),[1 ncunk]).*jb_re);
        tmp4=jb_im.'*(repmat(conj(y),[1 ncunk]).*jb_im);
    end
end

h2=4.0*real( tmp1+[ tmp2 tmp3.';tmp3 tmp4 ]  );

% Hessian of objective function
h_OBJ=h1+h2;

% HESSIAN OF THE SAR AND POWER CONSTRAINTS
if solver==1  % Matlab solver
    LagMult=lambda.ineqnonlin;
else
    LagMult=lambda;
end

% % Hessian of the constraints (SAR)
% subh=zeros(system.ncoils,system.ncoils);
% nconst=sar.nvop+1;
% subh=2.0*squeeze( sum( sar.sarmats .* repmat(LagMult(1:nconst),[1 nchannels nchannels]),1 ) );
% ch=kron( eye(spokes.nspokes),subh ).*(spokes.C);
% h_SAR=[real(ch) -imag(ch)';imag(ch)' real(ch)];

% Hessian of the constraints (max. power)
% tmp=LagMult(nconst+1:nconst+ncunk);
tmp=LagMult(1:ncunk);
tmp=repmat(tmp,[2 1]); % real and imaginary
h_MAX_POW=2.0*diag(tmp);

% h=h_OBJ + h_SAR + h_MAX_POW;
h=h_OBJ + h_MAX_POW;
% h=h_OBJ;





















