



function [f df]=f0_STA_SP(x_gr,adj,opt,spokes,system,do_mls,fac,SP_mean)
% function [f df]=f0_STA_SP(x_gr,adj,opt,spokes,system,do_mls,fac,SP_mean)

global A;
global spcoords;

gamma=42.576*1e6;

nchannels=size(adj.b1maps,1);
ncunk=spokes.nspokes*nchannels;

% retrieve RF
xrf=x_gr(2*spokes.nspokes+1:end);
rf=xrf(1:ncunk) + 1j*xrf(ncunk+1:end);

% compute system matrix if necessary
if sum( spcoords~=x_gr(1:2*spokes.nspokes,1)>0 )
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

% objective function
MT=A*rf;
if do_mls==1
    y=abs(MT) - abs(adj.target_STA);
else
    y=MT - adj.target_STA;
end
f=y'*y;

% gradient
if nargout>1    
    df=zeros(2*spokes.nspokes + 2*ncunk,1);
    
    % wrt spokes coordinates
    rf_diag=zeros( system.ncoils*spokes.nspokes,spokes.nspokes );
    for i=1:spokes.nspokes
        indrange=(i-1)*nchannels + 1 : i*nchannels;
        rf_diag( indrange,i )=rf(indrange,1);
    end
    tmp=2*pi*1j*A*rf_diag;  % Npix * Nspokes
    dM_dksx=tmp.*repmat( adj.nonzeropixels(:,4),[1 spokes.nspokes] ) / fac;
    dM_dksy=tmp.*repmat( adj.nonzeropixels(:,5),[1 spokes.nspokes] ) / fac;
    
    if do_mls==1
        a=2*real( dM_dksx'*( y.*exp(1j*angle(MT)) ) );
        b=2*real( dM_dksy'*( y.*exp(1j*angle(MT)) ) );
    else
        a=2*real( dM_dksx'*y );
        b=2*real( dM_dksy'*y ); 
    end
    df( 1:2:2*spokes.nspokes,1 )=a(end:-1:1,1);
    df( 2:2:2*spokes.nspokes,1 )=b(end:-1:1,1);
      
    % wrt RF
    if do_mls==1
        cdf=2*A'*( y.*exp(1j*angle(MT)) );
    else
        cdf=2.0*A'*y;        
    end
    df( end-2*ncunk+1:end,1 )=[ real(cdf);imag(cdf) ];
end


% force last spoke location @ (0,0)
df( (spokes.nspokes-1)*2+1 )=0;
df( (spokes.nspokes-1)*2+2 )=0;












