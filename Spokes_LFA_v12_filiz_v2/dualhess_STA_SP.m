


function h=dualhess_STA_SP(x_gr,lambda,fac,adj,spokes,sar,system,opt,SP_mean,do_mls,solver)
% function h=dualhess_STA_SP(x_gr,lambda,fac,adj,spokes,sar,system,opt,SP_mean,do_mls,solver)

global A;
global spcoords;

ncunk=spokes.nspokes * system.ncoils;

% retrieve RF
xrf=x_gr(2*spokes.nspokes+1:end);
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);

% compute system matrix if necessary
if spcoords~=x_gr(1:2*spokes.nspokes,1)
    spcoords=x_gr(1:2*spokes.nspokes,1);    
    % spokes locations
    for i=1:spokes.nspokes
        spokes.spcoords(i,1)=x_gr( (i-1)*2+1,1 )/fac;
        spokes.spcoords(i,2)=x_gr( (i-1)*2+2,1 )/fac;
    end    
    % system matrix
    spokes=prepare_pulse(spokes,opt,system);
    A=compute_system_matrix(adj,system,spokes,0,SP_mean);
end
MT=A*rf;


% HESSIAN OF OBJECTIVE (RF) -- 5%
if do_mls==1  
    y=abs(MT) - abs(adj.target_STA);
    z=sqrt( abs(adj.target_STA)./abs(MT) ) .* exp(1j*angle(MT));
    w=y./abs(MT);
    
    z( abs(MT)==0 )=0;
    w( abs(MT)==0 )=0;
    
    % X=A'*diag(z);
    X=A'.*repmat(z.',[ncunk 1]);  % faster    
    % Y=A'*diag(w)*A;
    Y=( A'.*repmat(w.',[ncunk 1]) )*A;  % faster
    
    h11=2*real(X)*real(X)' + 2*real(Y);
    h12=2*real(X)*imag(X)' - 2*imag(Y);
    h21=h12';
    h22=2*imag(X)*imag(X)' + 2*real(Y);
    
    h_RF=[ h11 h12;h21 h22 ];
else
    y=MT - adj.target_STA;
    ch=2.0*A'*A;
    h_RF=[real(ch) -imag(ch);
        imag(ch) real(ch)];
end


% HESSIAN OF OBJECTIVE (SPOKES COORDINATES) -- 10% seconds
rf_diag=zeros( system.ncoils*spokes.nspokes,spokes.nspokes );
for i=1:spokes.nspokes
    indrange=(i-1)*system.ncoils + 1 : i*system.ncoils;
    rf_diag( indrange,i )=rf(indrange,1);
end
tmp=A*rf_diag;  % Npix * Nspokes

dM_dkx=2*pi*1j*tmp.*repmat( adj.nonzeropixels(:,4),[1 spokes.nspokes] )/fac;
dM_dky=2*pi*1j*tmp.*repmat( adj.nonzeropixels(:,5),[1 spokes.nspokes] )/fac;

d2M_dkx_dkx=-4*pi^2*tmp.*repmat( adj.nonzeropixels(:,4).^2,[1 spokes.nspokes] )/fac^2;
d2M_dky_dky=-4*pi^2*tmp.*repmat( adj.nonzeropixels(:,5).^2,[1 spokes.nspokes] )/fac^2;
d2M_dkx_dky=-4*pi^2*tmp.*repmat( adj.nonzeropixels(:,4).*adj.nonzeropixels(:,5),[1 spokes.nspokes] )/fac^2;

if do_mls==1   
    z=abs(adj.target_STA)./abs(MT);
    z( abs(MT)==0 )=0;
    w=exp( 1j*angle(MT) );
    
    a=2*real( ( y.*w )'*d2M_dkx_dkx )'; a=a(end:-1:1);
    b=2*real( ( y.*w )'*d2M_dky_dky )'; b=b(end:-1:1);
    c=2*real( ( y.*w )'*d2M_dkx_dky )'; c=c(end:-1:1);

    z2=repmat( z,[1 spokes.nspokes] );    
    w2=repmat( w,[1 spokes.nspokes] );
    dyw_dkx=dM_dkx - z2.*( dM_dkx-w2.*real(conj(dM_dkx).*w2) );  % d_{ y*exp(1j*angle(MT)) }/dkx
    dyw_dky=dM_dky - z2.*( dM_dky-w2.*real(conj(dM_dky).*w2) );  % d_{ y*exp(1j*angle(MT)) }/dky
    
    d=2*real( dM_dkx'*dyw_dkx ); d=d(end:-1:1,end:-1:1);
    e=2*real( dM_dkx'*dyw_dky ); e=e(end:-1:1,end:-1:1);
    f=2*real( dM_dky'*dyw_dky ); f=f(end:-1:1,end:-1:1);
else
    a=2*real( y'*d2M_dkx_dkx )'; a=a(end:-1:1);
    b=2*real( y'*d2M_dky_dky )'; b=b(end:-1:1);
    c=2*real( y'*d2M_dkx_dky )'; c=c(end:-1:1);
    
    d=2*real( dM_dkx'*dM_dkx ); d=d(end:-1:1,end:-1:1);
    e=2*real( dM_dkx'*dM_dky ); e=e(end:-1:1,end:-1:1);
    f=2*real( dM_dky'*dM_dky ); f=f(end:-1:1,end:-1:1);
end
h_SC=zeros(2*spokes.nspokes);
h_SC( 1:2:end,1:2:end )=d + diag(a);  % dkx_dkx
h_SC( 2:2:end,1:2:end )=e + diag(c);  % dky_dkx
h_SC( 1:2:end,2:2:end )=e + diag(c);  % dkx_dky
h_SC( 2:2:end,2:2:end )=f + diag(b);  % dky_dky


% HESSIAN OF OBJECTIVE (SPOKES COORDINATES & RF) -- 8%
dA_dkx=2*pi*1j*A.*repmat( adj.nonzeropixels(:,4),[1 ncunk] )/fac;
dA_dky=2*pi*1j*A.*repmat( adj.nonzeropixels(:,5),[1 ncunk] )/fac;

if do_mls==1
    a=2*A'*dyw_dkx; a=a(:,end:-1:1); % Ncunk * Nspokes
    b=2*A'*dyw_dky; b=b(:,end:-1:1);    

    c=2*dA_dkx'*( y.*w ); % Ncunk * 1
    d=2*dA_dky'*( y.*w );
else
    a=2*A'*dM_dkx; a=a(:,end:-1:1); % Ncunk * Nspokes
    b=2*A'*dM_dky; b=b(:,end:-1:1);    
    
    c=2*dA_dkx'*y; % Ncunk * 1
    d=2*dA_dky'*y;
end
h_SC_RF=zeros( 2*ncunk,2*spokes.nspokes );
h_SC_RF( :,1:2:end )=[ real(a);imag(a) ];
h_SC_RF( :,2:2:end )=[ real(b);imag(b) ];

tmp=zeros( ncunk,2*spokes.nspokes );
for i=1:spokes.nspokes
    ind=(i-1)*system.ncoils+1:i*system.ncoils;
    tmp( ind,2*(i-1)+1:2*i )=[ d(ind) c(ind) ];  % the x,y order is reversed here so that the "tmp=tmp(:,end:-1:1);" operation 
                                                 % below reverses the spokes indices but not x and y.
end
tmp=tmp(:,end:-1:1);
h_SC_RF=h_SC_RF + [ real(tmp);imag(tmp) ];  % diagonal terms (when the spoke index of the RF=spoke index of the spoke coord derivative)


h_OBJ=[ h_SC h_SC_RF';
    h_SC_RF h_RF];


if solver==1  % Matlab solver
    LagMult=lambda.ineqnonlin;
else
    LagMult=lambda;
end

% % Hessian of the constraints (SAR)
% subh=zeros(system.ncoils,system.ncoils);
% nconst=sar.nvop+1;
% for n=1:nconst
%     subh(:,:)=subh(:,:) + 2.0*LagMult(n)*squeeze( sar.sarmats(n,:,:) );    
% end
% ch=kron( eye(spokes.nspokes),subh ).*(spokes.C);
% h=h+[real(ch) -imag(ch)';imag(ch)' real(ch)];

% Hessian of the constraints (max. power)
% tmp=LagMult(nconst+1:nconst+ncunk);
tmp=LagMult(1:ncunk);
tmp=[tmp;tmp]; % real and imaginary
h_CON=[ zeros( 2*spokes.nspokes,2*spokes.nspokes+2*ncunk );
    zeros( 2*ncunk,2*spokes.nspokes ) 2.0*diag(tmp) ];

h=h_OBJ + h_CON;
















